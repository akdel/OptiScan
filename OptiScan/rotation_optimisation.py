from scipy import ndimage
from scipy import signal
import numpy as np
from skimage.morphology import disk
from OptiScan.signal_match import Matcher
from OptiScan.align import normalized_correlation as ncorr
import matplotlib.pyplot as plt
from OptiScan.transformation import rotate as rotateit


def rotate(image, angle):
    return rotateit(image, np.deg2rad(angle))

# def rotate(image, angle):
#     return ndimage.rotate(image, angle, reshape=False)

def white_tophat_to_image(image_array, disk_radius=6):
    """
    Performs white tophat by using a disk structure. (Wrapper for ndimage.white_tophat + disk)
    :param image_array: Image numpy array; 2d numpy array.
    :param disk_radius: Radius for the disk structure; int.
    :return: White tophat result. 2D-Numpy array
    """
    struct = disk(disk_radius)
    return ndimage.white_tophat(image_array, structure=struct)


def get_xsum(image_array):
    """
    Sums the values in the array perpendicularly.
    :param image_array: Image numpy array; 2d numpy array.
    :return: 1D Numpy array.
    """
    return np.sum(image_array, axis=0)


def get_peaks_from_xsum(xsum, number=15):
    """
    Obtains peaks from bionano image xsum.
    :param xsum: Bionano image xsum; 1D-Numpy array.
    :return: Sorted tuples of (value, distance/index)
    """
    distances = signal.find_peaks_cwt(xsum, np.array([1]), min_snr=2)
    if distances is []:
        return None
    else:
        value_and_distance_pairs = [(xsum[x], x) for x in distances]
        return sorted(value_and_distance_pairs, reverse=True)[:number]


def get_average_peak_values(xsum):
    """
    :param xsum:
    :return:
    """
    peaks = get_peaks_from_xsum(xsum)
    if peaks:
        return np.mean([x[0] for x in peaks])
    else:
        return 0


def get_angle_range(_from, _to, space):
    """
    :param _from:
    :param _to:
    :param space:
    :return:
    """
    while _from <= _to:
        yield _from
        _from += space


def rotate_image_and_return_optimisation_value(image, rotation_angle, mid_xsum=False, single_peak=True):
    if not mid_xsum:
        image_xsum = ndimage.gaussian_filter1d(get_xsum(image), sigma=1)
        if abs(rotation_angle) <= 0.000000001:
            rotated_image_xsum = ndimage.grey_dilation(get_xsum(image), structure=np.ones((8)))
        else:
            rotated_image_xsum = ndimage.grey_dilation(get_xsum(rotate(image, rotation_angle)),
                                                       structure=np.ones((8)))
    else:
        im_shape = image.shape
        image_xsum = ndimage.gaussian_filter1d(
            get_xsum(image[:, (im_shape[1] / 2) - mid_xsum:(im_shape[1] / 2) + mid_xsum]),
            sigma=1)
        if abs(rotation_angle) <= 0.000000001:
            rotated_image_xsum = ndimage.grey_dilation(
                get_xsum(image[:, (im_shape[1] / 2) - mid_xsum:(im_shape[1] / 2) + mid_xsum]),
                structure=np.ones((8)))
        else:
            rotated_image_xsum = ndimage.grey_dilation(get_xsum(
                rotate(image[:, (im_shape[1] / 2) - mid_xsum:(im_shape[1] / 2) + mid_xsum], rotation_angle)), structure=np.ones((8)))
    if not single_peak:
        corr = Matcher(image_xsum, rotated_image_xsum)
        corr.convolve_signals_and_get_match_info()
        reference_peaks = get_peaks_from_xsum(corr.match[0])
        reference_peaks = [x for x in reference_peaks if x[0] >= 1000]
        if reference_peaks == []:
            return None
        original_image_peak_value = np.mean([x[0] for x in reference_peaks])
        reference_peaks = [x[1] for x in reference_peaks]
        rotation_image_peak_values =[corr.match[1][x] for x in reference_peaks]
        return np.mean(rotation_image_peak_values)/original_image_peak_value
    elif single_peak:
        return np.max(rotated_image_xsum)


def get_peak_averages_in_rotation_range(image, _from=-0.1, _to=0.1, space=0.05):
    """

    :param image:
    :param _from:
    :param _to:
    :param space:
    :return:
    """
    rotation_angles = [x for x in get_angle_range(_from, _to, space)]
    return [(rotate_image_and_return_optimisation_value(image, x, single_peak=True), x)
            for x in rotation_angles]


def get_optimal_rotation(image, _from=-0.1, _to=0.1, initial_space=0.1, final_space=0.01, saphyr=False):
    """
    :param image:
    :param _from:
    :param _to:
    :param initial_space:
    :param final_space:
    :return:
    """
    if saphyr:
        _from *= 10
        _to *= 10
        initial_space *= 2
    values_with_angles = get_peak_averages_in_rotation_range(image, _from=_from, _to=_to, space=initial_space)
    try:
        max_pair = max(values_with_angles)
    except TypeError:
        return [(1, 0)]
    rotation = max_pair[1]
    return get_peak_averages_in_rotation_range(image, _from=rotation-initial_space, _to=rotation+initial_space,
                                               space=final_space)


def rotate_with_optimal_rotation(image, _from=-0.1, _to=0.1, initial_space=0.05, final_space=0.01, saphyr=False):
    angle = max(get_optimal_rotation(image, _from=_from, _to=_to,
                                     initial_space=initial_space, final_space=final_space, saphyr=saphyr))[1]
    return rotate(image, angle), angle


def get_1d_bottom(image, saphyr=False):
    current_max = np.max(np.sum(image[-6:-1, :], axis=0))
    n = -1
    if saphyr:
        return np.mean(image[-200:], axis=0)
    else:
        while (current_max <= 600*5) and (n != -50):
            n -= 1
            current_max = np.max(np.sum(image[n-5:n, :], axis=0))
        return np.sum(image[n - 5:, :], axis=0)


def get_1d_top(image, saphyr=False):
    current_max = np.max(np.sum(image[0:5, :], axis=0))
    n = 0
    if saphyr:
        return np.mean(image[:200], axis=0)
    else:
        while (current_max <= 300*5) and (n != 50):
            n += 1
            current_max = np.max(np.sum(image[n:n+5, :], axis=0))
        return np.sum(image[n:n+5, :], axis=0)


def get_2d_bottom(image, saphyr=False):
    if saphyr:
        return image[-600:, :]
    return image[-120:, :]


def get_2d_top(image, saphyr=False):
    # print("Saphyr:",saphyr)
    if saphyr:
        return image[:600, :]
    return image[:120, :]


def x_shift_for_bottom_image(top_image, bottom_image, debug=False, saphyr=False):
    """
    :param top_image:
    :param bottom_image:
    :return:
    """
    bottom_image_top = get_1d_top(bottom_image, saphyr=saphyr)
    top_image_bottom = get_1d_bottom(top_image, saphyr=saphyr)
    if debug:
        plt.plot(top_image_bottom)
        plt.plot(bottom_image_top)
        plt.show()
        plt.plot(signal.correlate(top_image_bottom, bottom_image_top))
        plt.show()
    answer = ((np.argmax(signal.correlate(top_image_bottom, bottom_image_top)))
              - (-1 + len(bottom_image_top))) * -1
    if answer < 0:
        return ((np.argmax(signal.correlate(bottom_image_top, top_image_bottom)))
                - (-1 + len(bottom_image_top))) * -1, "top"
    else:
        return answer, "bottom"


def x_shift_image_while_keeping_default_xshape(shape, image, shift):
    if shift == 0:
        return image
    if image.shape[1] == shape[1]:
        new_frame = np.zeros(shape)
        new_frame[:, :-shift] = image[:, shift:]
        return new_frame
    else:
        return None


def x_shift_and_merge(top_image, bottom_image, shift_value, y_shift=False, return_y_shift=False, prey_shift=None, saphyr=False):
    # print(top_image.shape, bottom_image.shape)
    if top_image.shape[1] != bottom_image.shape[1]:
        return None
    if shift_value[1] == "bottom":
        bottom_image = x_shift_image_while_keeping_default_xshape(bottom_image.shape, bottom_image, shift_value[0])
    elif shift_value[1] == "top":
        top_image = x_shift_image_while_keeping_default_xshape(top_image.shape, top_image, shift_value[0])
    if y_shift:
        if prey_shift:
            if prey_shift is not 0:
                top_image = top_image[: -1 * prey_shift]
        else:
            top_bottom = get_2d_bottom(top_image, saphyr=saphyr)
            bottom_top = get_2d_top(bottom_image, saphyr=saphyr)
            try:
                _y = get_yshift(top_bottom, bottom_top, saphyr=saphyr)
            except ZeroDivisionError:
                _y = 0
            # _y = 0
            if _y != 0:
                top_image = top_image[: -_y]
    if return_y_shift:
        return np.concatenate((top_image, bottom_image)), _y
    return np.concatenate((top_image, bottom_image))


def x_shift_list_of_frames(list_of_frames_in_order, additional_set=None, y_shift=False, saphyr=False):
    current_frame = list_of_frames_in_order[0]
    if additional_set:
        current_additional_frame = additional_set[0]
    for i in range(1, len(list_of_frames_in_order), 1):
        shift_value = x_shift_for_bottom_image(current_frame, list_of_frames_in_order[i], saphyr=saphyr)
        current_frame, _y = x_shift_and_merge(current_frame, list_of_frames_in_order[i], shift_value, y_shift=y_shift,
                                              return_y_shift=True, saphyr=saphyr)
        if additional_set:
            current_additional_frame = x_shift_and_merge(current_additional_frame, additional_set[i], shift_value,
                                                         y_shift=True, prey_shift=_y)
    if additional_set:
        return current_frame, current_additional_frame
    else:
        return current_frame


def get_absolute_difference(array1, array2):
    """

    :param array1:
    :param array2:
    :return:
    """
    array1 = np.array(array1, dtype=int)
    array2 = np.array(array2, dtype=int)
    difference = array1 - array2
    absolute_difference = np.array([abs(x) for x in difference], dtype=int)
    return absolute_difference


def walk_up_and_get_absolute_differences(top_image_bottom, bottom_image_top):
    """

    :param top_image_bottom:
    :param bottom_image_top:
    :return:
    """
    bottom_1d_top = get_1d_top(bottom_image_top)
    for i in range(len(top_image_bottom))[::-1]:
        yield get_absolute_difference(bottom_1d_top, top_image_bottom[i, :])


def get_molecule_pairs(top_image_bottom, bottom_image_top):
    """
    :param top_image_bottom:
    :param bottom_image_top:
    :return:
    """
    top_bottom_peaks = [x[1] for x in get_peaks_from_xsum(get_xsum(top_image_bottom), number=30) if x[0] > 500]
    for xpeak in top_bottom_peaks:
        yield top_image_bottom[:, xpeak], bottom_image_top[:, xpeak]


def collapse_fields_into_2d_array(pairs):
    top_bottom_array = np.stack([x[0] for x in pairs])
    bottom_top_array = np.stack([x[1] for x in pairs])
    return top_bottom_array, bottom_top_array


def _vector_function_1(cell1, cell2):
    if cell1 == cell2:
        return 1
    else:
        return 0


def molecule_arrays_to_masks(top_bottom_array, bottom_top_array):
    top = np.where(top_bottom_array >= 300, 1, np.where(300 > top_bottom_array, 0, top_bottom_array))
    bot = np.where(bottom_top_array >= 300, 1, np.where(300 > bottom_top_array, 0, bottom_top_array))
    return top, bot


def vectorized_function_for_mask_difference(array1, array2):
    vectorized = np.vectorize(_vector_function_1)
    return vectorized(array1, array2)


def slide_and_score(array1, array2):
    for i in range(1, array1.shape[1], 1):
        array1_sliced = array1[:, -i:]
        array2_sliced = array2[:, :i]
        array1_masked, array2_masked = molecule_arrays_to_masks(array1_sliced, array2_sliced)
        overlap_array = vectorized_function_for_mask_difference(array1_masked, array2_masked)
        # plt.imshow(array2_masked, interpolation="nearest");plt.show()
        # plt.imshow(array1_masked, interpolation="nearest");plt.show()
        # plt.imshow(overlap_array, interpolation="nearest");plt.show()
        # plt.cla()
        yield len(np.where(overlap_array == 0)[0])**1.5/float(overlap_array.shape[0]*overlap_array.shape[1])


def top_bottom_to_slide_scores(top_image_bottom, bottom_image_top, return_score=False):
    pairs = [x for x in get_molecule_pairs(top_image_bottom, bottom_image_top)]
    xmed = (650*top_image_bottom.shape[0])/10
    pairs = [(x, y) for x, y in pairs if sum(x) >= xmed or sum(y) >= xmed]
    if len(pairs) == 0:
        return 0
    _2d_arrays = collapse_fields_into_2d_array(pairs)
    if return_score:
        return [x for x in slide_and_score(_2d_arrays[0], _2d_arrays[1])]
    else:
        return np.argmin([x for x in slide_and_score(_2d_arrays[0], _2d_arrays[1])])


def y_shift_by_mask_sliding(top_bottom_image, bottom_top_image):
    scores = top_bottom_to_slide_scores(top_bottom_image, bottom_top_image)
    if np.max(scores) < 0:
        return 0
    scores = ndimage.white_tophat(scores, structure=np.ones((6)))
    return np.argmax(scores)


def y_shift_for_list_of_frames(list_of_frames):
    for i in range(0, len(list_of_frames) - 1, 1):
        top_bottom_image = get_2d_bottom(list_of_frames[i])
        bottom_top_image = get_2d_top(list_of_frames[i+1])
        y_shift = y_shift_by_mask_sliding(top_bottom_image, bottom_top_image)
        yield y_shift


def merging_with_rotation_optimisation_and_xshift(list_of_frames, additional_set=None, y_shift=True, tophat=True,
                                                  magnification_optimisation=True, saphyr=False):
    """
    This is the pipeline function and is used for bionano column alignment..
    :param list_of_frames:
    :param additional_set:
    :param y_shift:
    :param tophat:
    :param magnification_optimisation:
    :return:
    """
    if tophat:
        list_of_frames = [ndimage.white_tophat(x, structure=disk(7)) for x in list_of_frames] # 6 for irys
        if additional_set:
            additional_set = [ndimage.white_tophat(x, structure=disk(9)) for x in additional_set]
    list_of_frames_with_angles = [rotate_with_optimal_rotation(x, saphyr=saphyr) for x in list_of_frames]
    list_of_frames = [x[0].astype(float) for x in list_of_frames_with_angles]
    angles = [x[1] for x in list_of_frames_with_angles]
    if additional_set:
        additional_input = [rotate(additional_set[i], angles[i]).astype(float)
                            for i in range(len(additional_set))]
        if magnification_optimisation:
            additional_mag_input = [get_optimal_zoom_and_obtain_new_image(additional_input[i], list_of_frames[i]) for i
                                    in range(len(additional_input))]
            return x_shift_list_of_frames(list_of_frames, additional_set=additional_mag_input, y_shift=y_shift, saphyr=saphyr)
        return x_shift_list_of_frames(list_of_frames, additional_set=additional_input, y_shift=y_shift, saphyr=saphyr)
    else:
        return x_shift_list_of_frames(list_of_frames, additional_set=additional_set, y_shift=y_shift, saphyr=saphyr)


def get_yshift2(top_image_bottom, bottom_image_top, return_score=False):
    pairs = [(top_image_bottom[:,i], bottom_image_top[:,i]) for i in range(0, top_image_bottom.shape[1], 1)]
    if len(pairs) == 0:
        return 0
    xmed = np.median([sum(x) for x, y in pairs])
    ymed = np.median([sum(y) for x, y in pairs])
    filtered_pairs = [(x,y) for x, y in pairs if (sum(x)/1.5 >= xmed) or (sum(y)/1.5 >= ymed)]
    corrs = [np.correlate(y, x, mode="full") for x, y in filtered_pairs]
    if not corrs:
        if return_score:
            return 0,0
        else:
            return 0
    corr_sum = sum(corrs)
    if return_score:
        return corr_sum, corrs
    return np.argmax(corr_sum[:60])


def get_yshift(top_image_bottom, bottom_image_top, debug=True, saphyr=False):
    pairs = [(top_image_bottom[:,i], bottom_image_top[:,i]) for i in range(0, top_image_bottom.shape[1], 1)]
    xmed = np.median([sum(x) for x, y in pairs])
    ymed = np.median([sum(y) for x, y in pairs])
    if saphyr:
        filtered_pairs = [(x,y) for x, y in pairs if (sum(x) >= xmed * 4) or (sum(y) >= ymed * 4)]
        corrs = np.array([ncorr(y, x, limit=25) for x, y in filtered_pairs], dtype=float)
    else:
        filtered_pairs = [(x,y) for x, y in pairs if (sum(x) >= xmed * 4) or (sum(y) >= ymed * 4)]
        corrs = np.array([ncorr(y, x, limit=8) for x, y in filtered_pairs], dtype=float) #limit=12
    if len(filtered_pairs) == 0:
        return 0
    corr_sum = np.sum(corrs, axis=0)
    if debug:
        plt.imshow(top_image_bottom)
        plt.show()
        plt.imshow(bottom_image_top)
        plt.show()
        [plt.plot(corr) for corr in corrs]
        plt.show()

        plt.plot(corr_sum)
        plt.show()
    else:
        pass
    if saphyr:
        return np.argmax(corr_sum[:300])
    else:
        return np.argmax(corr_sum[:40]) #limit=60



def zoom_out_and_center_on_original(image, zoom_out_ratio, shift):
    if zoom_out_ratio == 1:
        return image
    original_size = image.shape[0]
    shrank = ndimage.zoom(image, zoom_out_ratio)
    new_shrank = np.zeros(shrank.shape, dtype=shrank.dtype)
    if shift == 0:
        new_shrank = shrank
    elif shift < 0:
        new_shrank[:,:shift] = shrank[:,-1*shift:]
    else:
        new_shrank[:,shift:] = shrank[:,:-shift] # fix
    shrank = new_shrank
    shrank_size = shrank.shape[0]
    size_diff = original_size - shrank_size
    pad_width = size_diff / 2
    if size_diff % 2 != 0:
        shrank = np.concatenate((shrank, np.zeros((shrank.shape[0], 1))), axis=1)
        shrank = np.concatenate((shrank, np.zeros((1, shrank.shape[1]))), axis=0)
    if zoom_out_ratio > 1:
        print("Warning! incorrect ratio, it should be less than 1")
        return image
    else:
        return np.pad(shrank, int(pad_width), "constant")


def get_corr_score_for_zoom(image_xsum, ref_image_xsum, zoom_out_ratio):
    image_xsum = ndimage.zoom(image_xsum, zoom_out_ratio)
    corr = signal.correlate(ref_image_xsum, image_xsum)
    # print(image_xsum.shape, ref_image_xsum.shape)
    # plt.plot(corr)
    # plt.show()
    shift_idx = np.argmax(corr)
    max_corr = corr[shift_idx]
    # print(max_corr, shift_idx)
    return max_corr, shift_idx - image_xsum.shape[0]
    


def get_optimal_magnification_for_overlay(image, ref_image, _start=0.990, _to=0.999, _space=0.001):
    zoom_values = get_angle_range(_start, _to, _space)
    image_xsum = get_xsum(image)
    ref_image_xsum = get_xsum(ref_image)
    mags = [(get_corr_score_for_zoom(image_xsum, ref_image_xsum, x), x) for x in zoom_values]
    max_mag =  max(mags)
    return max_mag[1], max_mag[0][1]


def get_optimal_zoom_and_obtain_new_image(image, ref_image, _start=0.990, _to=0.999, _space=0.001):
    optimal_mag, shift = get_optimal_magnification_for_overlay(image, ref_image, _start=_start, _to=_to, _space=_space)
    return zoom_out_and_center_on_original(image, optimal_mag, shift)


def overlay_saphyr_columns(mol_col, label_col, _start=0.990, _to=1.01, _space=0.001):
    optimized_mol_col, angle =rotate_with_optimal_rotation(mol_col)
    rotated_label_col = rotate(label_col, angle)
    optimized_label_col = get_optimal_zoom_and_obtain_new_image(rotated_label_col, optimized_mol_col, _start=_start, _to=_to, _space=0.001)
    return optimized_mol_col, optimized_label_col

"""

def _correlate_molecule_pairs(molecule_pairs):
    for pair in molecule_pairs:
        corralation = signal.correlate(pair[0], pair[1])
        cor_le
        n = len(corralation)
        mol_len = cor_len/2
        print cor_len
         mol_len
        init_adjusted_corralation = [corralation[x]/float(x) for x in range(cor_len) if x <= mol_len]
        end_adjusted_corralation = [corralation[x]/float(mol_len-abs(mol_len-x)) for x in range(cor_len) if x > mol_len]
        yield init_adjusted_corralation + end_adjusted_corralation


def _corralate_molecules_from_top_and_bottom(top_image_bottom, bottom_image_top):
    molecule_pairs = [x for x in get_molecule_pairs(top_image_bottom, bottom_image_top)]
    return _correlate_molecule_pairs(molecule_pairs)


def _get_max_correlation_values(molecule_corralations):
    for corralation in molecule_corralations:
        yield int(np.argmax(corralation))


def get_max_corralation_median_of_molecule_pairs_from_top_and_bottom(top_image_bottom, bottom_image_top):
    molecule_corralations = _corralate_molecules_from_top_and_bottom(top_image_bottom, bottom_image_top)
    return np.median([x for x in _get_max_correlation_values(molecule_corralations)])


def get_optimal_start_points_for_y_shift_estimation(top_image_bottom, bottom_image_top):
    list_of_abs_differences = [x for x in walk_up_and_get_absolute_differences(top_image_bottom, bottom_image_top)]
    minimum_difference_y_shift = np.argmin([sum(x) for x in list_of_abs_differences])
    return minimum_difference_y_shift

def extend_bottom_top_until_high_difference(y_shift, top_image_bottom, bottom_image_top):
    inverted_top_image_bottom = top_image_bottom[::-1]
    bottom_window = inverted_top_image_bottom[:y_shift]
    return get_optimal_start_points_for_y_shift_estimation(bottom_image_top[::-1], bottom_window), y_shift

"""

