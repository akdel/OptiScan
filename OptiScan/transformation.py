import numpy as np
import numba as nb

canvas_dtype = np.dtype({"names": ["red", "green", "blue"], "formats": ["f4", "f4", "f4"]})
ppm_dtype = np.dtype({"names": ["red", "green", "blue"], "formats": ["uint8", "uint8", "uint8"]})
tuple_dtype = np.dtype({"names": ["x", "y", "z", "w"], "formats": ["f8", "f8", "f8", "f8"]})


@nb.njit
def pad_for_ppm(val):
    if val <= 0:
        return int(0)
    elif val >= 1:
        return int(255)
    else:
        return int(256 * val)


@nb.njit
def convert_canvas_to_ppm(canvas):
    ppm = np.zeros(canvas.shape, dtype=ppm_dtype)
    for i in nb.prange(ppm.shape[0]):
        for j in range(ppm.shape[1]):
            ppm.red[i,j] = pad_for_ppm(canvas.red[i,j])
            ppm.green[i, j] = pad_for_ppm(canvas.green[i, j])
            ppm.blue[i, j] = pad_for_ppm(canvas.blue[i, j])
    return ppm


@nb.njit
def projectile(vector, velocity, iter=100, gravity=1.0, wind=0.5):
    points = np.zeros((iter + 1, 4), dtype=np.float64)
    points[0] = vector
    for i in range(1, points.shape[0]):
        points[i] = velocity + points[i-1]
        velocity[1] -= gravity
        velocity[0] += wind
        print(velocity)
    return points


@nb.njit
def points_within_envelope(points, xmax, ymax):
    for i in range(points.shape[0]):
        if points[i][0] > xmax:
            points[i][0] = xmax
        elif points[i][0] < 0:
            points[i][0] = 0
        if points[i][1] > ymax:
            points[i][1] = ymax
        elif points[i][1] < 0:
            points[i][1] = 0
    return points


@nb.njit
def draw_from_points(points, canvas):
    for i in nb.prange(points.shape[0]):
        x, y = points[i][0], points[i][1]
        canvas[x,y].red = 1
        canvas[x,y].blue = 1
        canvas[x,y].green = 1
    return canvas


def write_ppm_to_file(ppm, filename="./ppm.ppm"):
    f = open(filename, "w")
    header = "P3\n%s %s\n255\n" % ppm.shape
    f.write(header)
    rest = list()
    _ = [rest.append(" ".join(map(str,val))) for val in ppm.flatten()]
    f.write(" ".join(rest))
    f.write("\n")
    f.close()


@nb.njit
def create_translation_matrix(x, y, z):
    translation_matrix = np.identity(4)
    translation_matrix[:3, -1] = (x, y, z)
    return translation_matrix


@nb.njit
def create_scale_matrix(x, y, z):
    scale_matrix = np.identity(4)
    scale_matrix[0, 0] = x
    scale_matrix[1, 1] = y
    scale_matrix[2, 2] = z
    return scale_matrix


@nb.njit
def create_x_rotation_matrix(angle):
    matrix = np.identity(4)
    matrix[1, 1] = np.cos(angle)
    matrix[2, 2] = np.cos(angle)
    matrix[1, 2] = -np.sin(angle)
    matrix[2, 1] = np.sin(angle)
    return matrix


@nb.njit
def create_z_rotation_matrix(angle):
    matrix = np.identity(4)
    matrix[0, 0] = np.cos(angle)
    matrix[0, 1] = -np.sin(angle)
    matrix[1, 0] = np.sin(angle)
    matrix[1, 1] = np.cos(angle)
    return matrix


@nb.njit
def create_y_rotation_matrix(angle):
    matrix = np.identity(4)
    matrix[0, 0] = np.cos(angle)
    matrix[0, 2] = np.sin(angle)
    matrix[2, 0] = -np.sin(angle)
    matrix[2, 2] = np.cos(angle)
    return matrix


@nb.njit
def create_shear_matrix(xy, xz, yx, yz, zx, zy):
    matrix = np.identity(4)
    matrix[0, 1] = xy
    matrix[0, 2] = xz
    matrix[1, 0] = yx
    matrix[1, 2] = yz
    matrix[2, 0] = zx
    matrix[2, 1] = zy
    return matrix


class Transformation:
    def __init__(self, tuples):
        self.tuples = tuples.T

    def translate(self, x, y, z, replace=False):
        translation_matrix = create_translation_matrix(x, y, z)
        if replace:
            self.tuples = translation_matrix @ self.tuples
        else:
            return translation_matrix @ self.tuples

    def scale(self, x, y, z, replace=False):
        scale_matrix = create_scale_matrix(x, y, z)
        if replace:
            self.tuples = scale_matrix @ self.tuples
        else:
            return scale_matrix @ self.tuples

    def custom_transformation(self, matrix, replace=False):
        if replace:
            self.tuples = matrix @ self.tuples
        else:
            return matrix @ self.tuples

    def reflect_on_axis(self, axis_name, replace=False):
        try:
            assert "x" in axis_name or "y" in axis_name or "z" in axis_name
        except AssertionError:
            raise AssertionError("give 'x', 'y' or 'z' as axis name")
        if axis_name == "x":
            return self.scale(-1, 1, 1, replace=replace)
        if axis_name == "y":
            return self.scale(1, -1, 1, replace=replace)
        else:
            return self.scale(1, 1, -1, replace=replace)

    def rotate(self, axis_name, angle, replace=False):
        try:
            assert "x" in axis_name or "y" in axis_name or "z" in axis_name
        except AssertionError:
            raise AssertionError("give 'x', 'y' or 'z' as axis name")
        if axis_name == "x":
            r_matrix = create_x_rotation_matrix(angle)
        elif axis_name == "y":
            r_matrix = create_y_rotation_matrix(angle)
        else:
            r_matrix = create_z_rotation_matrix(angle)
        if replace:
            self.tuples = r_matrix @ self.tuples
        else:
            return r_matrix @ self.tuples

    def shear(self, xy, xz, yx, yz, zx, zy, replace=False):
        if replace:
            self.tuples = create_shear_matrix(xy, xz, yx, yz, zx, zy) @ self.tuples
        else:
            return create_shear_matrix(xy, xz, yx, yz, zx, zy) @ self.tuples


@nb.njit
def pixels_to_vectors(pixels):
    res = np.zeros((pixels.shape[0] * pixels.shape[1], 4), dtype=np.float64)
    for i in range(pixels.shape[0]):
        for j in range(pixels.shape[1]):
            res[(i * pixels.shape[1]) + j] = (j, i, pixels[i, j], 1)
    return res


@nb.njit
def vectors_to_pixels(vectors):
    length = int(np.sqrt(vectors.shape[0]))
    res = np.zeros((length, length), dtype=np.int64)
    for i in range(vectors.shape[0]):
        v = vectors[i]
        res[v[1], v[0]] = v[2]
    return res


@nb.jit(nopython=True, parallel=True)
def rotate(image, angle):
    x_len = image.shape[1]
    y_len = image.shape[0]
    tuples = pixels_to_vectors(image).T
    t_mat = create_translation_matrix(x_len//2, y_len//2, 0) @ \
            create_z_rotation_matrix(angle) @ \
            create_translation_matrix(-x_len//2, -y_len//2, 0)
    rotated = t_mat @ tuples
    transformed = points_within_envelope(rotated.T, x_len-1, y_len-1).astype(np.int64)
    return vectors_to_pixels(transformed)

if __name__ == "__main__":
    canvas = np.zeros((500,500), dtype=canvas_dtype)
    canvas[:,:] = -1
    points = projectile(np.array([0, 250, 0, 0], dtype=float), np.array([2, 5, 0, 1], dtype=float), iter=130, gravity=0.1, wind=0.025)
    points = points_within_envelope(points, 499, 499)
    t = Transformation(points)
    points = t.rotate("z", np.pi/10, replace=False)
    draw_from_points(points.astype(int), canvas)
    ppm = convert_canvas_to_ppm(canvas)
    write_ppm_to_file(ppm, "ppm4.ppm")
