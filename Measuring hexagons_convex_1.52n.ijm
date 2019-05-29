// BIC_Measuring Hexagons
// =======================================================
// Copyright 2018 BioImaging Center, University of Konstanz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version. 
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. If data produced 
// with this program or a derivative of this program is used for 
// publications the original authors should be acknowledged appropriately. 
//
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// =======================================================
//
//	Uses "FeatureJ Laplacian" written by Erik Meijering 
//	<https://imagej.net/ImageScience>
//
// =======================================================
//
// BioImaging Center, University of Konstanz <bioimaging@uni-konstanz.de>
//   Martin Stöckl <martin.stoeckl@uni-konstanz.de>
//
// =======================================================

/* returns first index of element in an array: returns -1 if element not found */
function findIndex(array, element) {
	for (i=0; i < lengthOf(array); i++) {
		if (array[i] == element) {
			return i;
		}
	}
	return -1;
}

/* distance between two points, pythagorean theorem */
function calculateDistance(x0, y0, x1, y1) {
	return sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));
}

/* helper function, returns an array of angles between an array of points and a common vertex */
function calculateAngles(x_corners, y_corners, x_vertex, y_vertex) {
	angles =  newArray(x_corners.length);
	for (c = 0; c < x_corners.length; c++){
		angle = calculateAngle(x_corners[c], y_corners[c], x_corners[(c + 1) % x_corners.length], y_corners[(c + 1) % x_corners.length], x_vertex, y_vertex);		
		angles[c] = angle;
	}
	return angles;	
}

/* calculates the angle between two points and a vertex point, law of cosines */
function calculateAngle(x0, y0, x1, y1, x_vertex, y_vertex) {
	length0 = calculateDistance(x0, y0, x_vertex, y_vertex);
	length1 = calculateDistance(x1, y1, x_vertex, y_vertex);
	length_opposite = calculateDistance(x0, y0, x1, y1);
	angle_cos = (pow(length0, 2) + pow(length1, 2) - pow(length_opposite, 2)) / (2 * length0 * length1);
	
	return acos(angle_cos) * 180 / PI;
}

/* sum of intensities along a profile, returns an array of sums for given angles */
// title --- image title
// x_center, y_center --- origin point of profile
// radius --- length of profile
// start_angle, stop_angle --- angle values (in degrees) between which the profiles are analysed
// normalize --- return real or normalized values (between 0 and 1)
function getProfilesAsArray(title, x_center, y_center, radius, start_angle, stop_angle, normalize) {
	
	profile_sums = newArray();
	if (stop_angle < start_angle) {
		stop_angle += 360;		
	}
	
	selectImage(title);		
	for (j = start_angle; j < stop_angle; j++) {
		x = x_center + radius * cos(j * PI / 180);
		y = y_center + radius * sin(j * PI / 180);
		
		makeLine(x_center, y_center, x, y);
		profile = getProfile();
				
		sum_profile = 0;
		for (n = 0; n < lengthOf(profile); n++) {
			sum_profile += profile[n];			
		}		
		profile_sums = Array.concat(profile_sums, sum_profile);		
	}
	/* to increase the response average number of angles given by SAMPLING_WINDOW */
	resampled_sums = newArray();
	for (slice_start = 0; slice_start <= (profile_sums.length - SAMPLING_WINDOW); slice_start += SAMPLING_WINDOW) {
		current_slice = Array.slice(profile_sums, slice_start, slice_start + SAMPLING_WINDOW);
		Array.getStatistics(current_slice, VOID, VOID, slice_mean, VOID);	
		resampled_sums = Array.concat(resampled_sums, slice_mean);
	}
	
	if (!normalize) {
		return resampled_sums;	
	}
	else {
		normalized_sums = newArray(resampled_sums.length);
		Array.getStatistics(resampled_sums, min_sum, max_sum, VOID, VOID);
		for (s = 0; s < resampled_sums.length; s++) {
			normalized_sums[s] = (resampled_sums[s] - min_sum) / (max_sum - min_sum);
		}		
		return normalized_sums;
	}
}

/* calculate edge distances from center, returns an array of distances */
// title --- image title
// x_center, y_center --- origin point of profile
// radius --- length of profile
// normalize --- return real or normalized values (between 0 and 1)
function getEdgesAsArray(title, x_center, y_center, radius, start_angle, stop_angle, normalize) {

	profile_edges = newArray();
	if (stop_angle < start_angle) {
		stop_angle += 360;		
	}
	
	selectImage(title);		
	for (j = start_angle; j < stop_angle; j++) {
		
		/* create end point outside the hexagon */
		x = x_center + radius * cos(j * PI / 180);
		y = y_center + radius * sin(j * PI / 180);	
			
		/* get coordinate of edge pixel and calculate distance from center */
		edge_hit = refineCornerDistance(x_center, y_center, x, y);		
		edge_distance = calculateDistance(x_center, y_center, edge_hit[0], edge_hit[1]);
		profile_edges = Array.concat(profile_edges, edge_distance);		
	}
	
	/* to increase the response average number of angles given by SAMPLING_WINDOW */
	resampled_sums = newArray();
	for (slice_start = 0; slice_start <= (profile_edges.length - SAMPLING_WINDOW); slice_start += SAMPLING_WINDOW) {
		current_slice = Array.slice(profile_edges, slice_start, slice_start + SAMPLING_WINDOW);
		Array.getStatistics(current_slice, VOID, VOID, slice_mean, VOID);	
		resampled_sums = Array.concat(resampled_sums, slice_mean);
	}

	if (!normalize) {
		return resampled_sums;	
	}
	else {
		normalized_sums = newArray(resampled_sums.length);
		Array.getStatistics(resampled_sums, min_sum, max_sum, VOID, VOID);
		for (s = 0; s < resampled_sums.length; s++) {
			normalized_sums[s] = (resampled_sums[s] - min_sum) / (max_sum - min_sum);
		}		
		return normalized_sums;
	}
} 

/* returns local maxima present in an array_int */
// find_central --- if true: return only one maximum, exclude maxima at array boundaries, returns -1 if no central maximum
function detectMaxima(array_int, find_central) {

	hits = newArray();
	current_value = array_int[0];
	
	/* determine a possible max at index 0 */
	if (array_int[0] < array_int[1]) {
		goes_up = true;
	}
	else {
		goes_up = false;
		hits = Array.concat(hits, current_value);
	}
	
	/* after a max, values start to decrease */
	for (i = 1; i < lengthOf(array_int); i++) {
		if (array_int[i] < current_value && goes_up){
			goes_up = false;
			hits = Array.concat(hits, current_value);			
		}
		if (array_int[i] >= current_value) {
			goes_up = true;
		}
		
		current_value = array_int[i];	
	}
	/* determine a possible max at the last index */
	if (goes_up) {
		hits = Array.concat(hits, current_value);
	}

	/* find the indices for the max intensities */
	indices_found = newArray(lengthOf(hits));
	last_index = -1;
	for (h = 0; h < lengthOf(hits); h++) {
		current_slice = Array.slice(array_int, last_index + 1);
		indices_found[h] = findIndex(current_slice, hits[h]) + last_index + 1;
		last_index = indices_found[h];
	}

	/* calculate distances between maxima */
	hit_distances = newArray(lengthOf(hits));	
		 
	for (h = 0; h < lengthOf(hits) - 1; h++) {
		hit_distances[h] = indices_found[h + 1] - indices_found[h];		
	}
	hit_distances[lengthOf(hits) - 1] = array_int.length - indices_found[lengthOf(hits) - 1] + indices_found[0];

	/* For debugging: Print identified maxima to log */	
	//print("Detected maxima:");
	//for (t = 0; t < lengthOf(hits); t++) {
	//	print("Int " + hits[t] + " Index " + indices_found[t], " Distance " + hit_distances[t]);
	//}
	

	/* return only the highest local max */
	// maxima closer than IGNORE_LIMIT to the array boundaries are ignored
	if (find_central) {
				
		central_maximum = newArray();
		central_max_int = newArray();
		
		for (m = 0; m < indices_found.length; m++) {
			if (indices_found[m] > IGNORE_LIMIT && ((array_int.length - indices_found[m]) > IGNORE_LIMIT)) {
				central_maximum = Array.concat(central_maximum, indices_found[m]);
				central_max_int = Array.concat(central_max_int, hits[m]);
			}
		}
		
		if (central_maximum.length > 0) {
			Array.getStatistics(central_max_int, min, max_int, mean, stdDev);
			max_int_index = findIndex(central_max_int, max_int);
			return central_maximum[max_int_index];
		}
		else {
			return -1;
		}
	}	
	
	/* identify close maxima with a separation smaller than MAX_HITDISTANCE and group them into one (indices are averaged) */
	separated_indices = newArray();
	separated_intensities = newArray();
	grouped_indices = newArray();
	grouped_intensities = newArray();

	/* if max distance smaller MAX_HITDISTANCE add to group, else add to and average current group, add to separated max, initialize new group, if no current group add directly to separated max */
	for (h = 0; h < lengthOf(hits); h++) {		
		if (hit_distances[h] <= MAX_HITDISTANCE) {
			grouped_indices = Array.concat(grouped_indices, indices_found[h]);
			grouped_intensities = Array.concat(grouped_intensities, hits[h]);
		}
		else if (hit_distances[h] > MAX_HITDISTANCE && lengthOf(grouped_indices) > 0) {
			grouped_indices = Array.concat(grouped_indices, indices_found[h]);
			grouped_intensities = Array.concat(grouped_intensities, hits[h]);
			Array.getStatistics(grouped_indices, VOID, VOID, mean_index, VOID);
			separated_indices = Array.concat(separated_indices, mean_index);
			Array.getStatistics(grouped_intensities, VOID, VOID, mean_int, VOID);
			separated_intensities = Array.concat(separated_intensities, mean_int);
			grouped_indices = newArray();
			grouped_intensities = newArray();				
		}
		else {
			separated_indices = Array.concat(separated_indices, indices_found[h]);
			separated_intensities = Array.concat(separated_intensities, hits[h]);
		}
	}

	/* clean up grouped max, check if max across array boundaries have to be grouped, add to separated max */
	if (lengthOf(grouped_indices) > 0) {
		if (hit_distances[h - 1] > MAX_HITDISTANCE) {
			Array.getStatistics(grouped_indices, VOID, VOID, mean_index, VOID);
			separated_indices = Array.concat(separated_indices, mean_index);
			Array.getStatistics(grouped_intensities, VOID, VOID, mean_int, VOID);
			separated_intensities = Array.concat(separated_intensities, mean_int);
		}
		else {
			if (lengthOf(separated_indices) > 0) {
				grouped_indices = Array.concat(grouped_indices, separated_indices[0] + 360 / SAMPLING_WINDOW);
				grouped_intensities = Array.concat(grouped_intensities, separated_intensities[0]);
				Array.getStatistics(grouped_indices, VOID, VOID, mean_index, VOID);
				if (mean_index >= 360 / SAMPLING_WINDOW) {
					mean_index -= 360 / SAMPLING_WINDOW;
				}
				separated_indices = Array.concat(separated_indices, mean_index);
				Array.getStatistics(grouped_intensities, VOID, VOID, mean_int, VOID);
				separated_intensities = Array.concat(separated_intensities, mean_int);
				separated_indices = Array.slice(separated_indices, 1);
				separated_intensities = Array.slice(separated_intensities, 1);			
			}
			else {
				Array.getStatistics(grouped_indices, VOID, VOID, mean_index, VOID);
				separated_indices = Array.concat(separated_indices, mean_index);
				Array.getStatistics(grouped_intensities, VOID, VOID, mean_index, VOID);
				separated_intensities = Array.concat(separated_intensities, mean_index);
			}
		}
	}
	
	/* For debugging: Print identified maxima to log */
	//print("After separation:");
	//for (t = 0; t < lengthOf(separated_indices); t++) {
	//	print("Int " + separated_intensities[t] + " Index " + separated_indices[t]);
	//}
		
	return separated_indices;		
}

/* for each corner angle, check if an opposite corner angle is there, if not add one */
// angles --- array of corner angles
function probeOppositeAngles(angles) {
		
	opp_angles = newArray();
	
	for (a = 0; a < angles.length; a++) {
		if (angles[a] >= 180 / SAMPLING_WINDOW) {
			opp_angle = angles[a] - 180 / SAMPLING_WINDOW;
		}
		else {
			opp_angle = angles[a] + 180 / SAMPLING_WINDOW;
		}
		
		/* check if there is already an angle closer than MIN_OPP_ANGLE to the new opposite angle */		
		angle_close = false;
		for (av = 0; av < angles.length; av++) {
			if (abs(angles[av] - opp_angle) < MIN_OPP_ANGLE / SAMPLING_WINDOW) {
				angle_close = true;
			}
		}
		if (!angle_close) {
			opp_angles = Array.concat(opp_angles, opp_angle);		
		}		
	}
	angles = Array.concat(angles, opp_angles);
	Array.sort(angles);
	return angles;	
}

/* refine edge distance from center, for a single point */
// x_center, y_center --- origin point of profile
// x_initial, y_initial --- initial point coordinates
function refineCornerDistance(x_center, y_center, x_initial, y_initial) {
	
	/* get angle from center, extend search range by 7 pixel */
	angle = atan2(y_initial - y_center, x_initial - x_center);
	x_ext = cos(angle) * 7;
	y_ext = sin(angle) * 7;	
	x_extended = x_initial + x_ext;
	y_extended = y_initial + y_ext;	
	
	/* identify step in intesities along a profile, by discrete differentitaion in sample window given by DIFFERENTIATION_WINDOW */
	selectImage("Working");		
	makeLine(x_center, y_center, x_extended, y_extended);	
	profile = getProfile();
	x_points = Array.getSequence(DIFFERENTIATION_WINDOW);
	limit = floor(DIFFERENTIATION_WINDOW / 2);
	diffed_profile = newArray(profile.length);					
	for (w = limit; w < profile.length - limit; w++) {					
		temp_array = Array.slice(profile, w - limit, w + limit + 1);	
		Fit.doFit("Straight Line", x_points, temp_array);
		diffed_profile[w] = Fit.p(1);		
	}
	/* identify index with maximum slope */
	Array.getStatistics(diffed_profile, min, max_differential, mean, stdDev);
	corner_dist = findIndex(diffed_profile, max_differential);

	/* calculate refined coordinates */
	coord_mod = newArray(2);
	coord_mod[0] = cos(angle) * corner_dist + x_center;
	coord_mod[1] = sin(angle) * corner_dist + y_center;
	
	return coord_mod;
}

/*  calculates angles between adjacent points from an sorted array of corner points, eliminate corners with an angle above MAX_ALLOWED_ANGLE */
// corners --- array of corner coordinates [[x_coords] + [y_coords]]
function removeCornersOnEdgesAngleLimit(corners) {

	x_corners_new = newArray();
	y_corners_new = newArray();
	x_corners = Array.slice(corners, 0, corners.length / 2);
	y_corners = Array.slice(corners, corners.length / 2);
	
	for (c = 0; c < x_corners.length; c++) {
		x0 = x_corners[c];
		y0 = y_corners[c];
		x_v = x_corners[(c + 1) % x_corners.length];
		y_v = y_corners[(c + 1) % y_corners.length];
		x1 = x_corners[(c + 2) % x_corners.length];
		y1 = y_corners[(c + 2) % y_corners.length];
		
		corner_angle = calculateAngle(x0, y0, x1, y1, x_v, y_v);
		if (corner_angle <= MAX_ALLOWED_ANGLE) {
			x_corners_new = Array.concat(x_corners_new, x_v);
			y_corners_new = Array.concat(y_corners_new, y_v);			
		}
	}
	corners_new = Array.concat(x_corners_new, y_corners_new);
	return corners_new;
}

/*  calculates angles between adjacent points from an sorted array of corner points, return the six corners with smallest angles */
// if not at least seven corners are given, fall back to removeCornersOnEdgesAngleLimit
// corners --- array of corner coordinates [[x_coords] + [y_coords]]
function removeCornersOnEdges(corners) {

	if (corners.length <= 12) {
		corners_new = removeCornersOnEdgesAngleLimit(corners);
	}
	else {
		x_corners_new = newArray(6);
		y_corners_new = newArray(6);
		x_corners = Array.slice(corners, 0, corners.length / 2);
		y_corners = Array.slice(corners, corners.length / 2);
		angles = newArray();
		corner_indices = newArray(6);
		
		for (c = 0; c < x_corners.length; c++) {
			x0 = x_corners[c];
			y0 = y_corners[c];
			x_v = x_corners[(c + 1) % x_corners.length];
			y_v = y_corners[(c + 1) % y_corners.length];
			x1 = x_corners[(c + 2) % x_corners.length];
			y1 = y_corners[(c + 2) % y_corners.length];
			
			corner_angle = calculateAngle(x0, y0, x1, y1, x_v, y_v);
			angles = Array.concat(angles, corner_angle);		
		}
		sorted_angles = Array.copy(angles);
		Array.sort(sorted_angles);	
		for (i = 0; i < 6; i++) {
			corner_indices[i] = (findIndex(angles, sorted_angles[i]) + 1) % x_corners.length;		
		}	
		Array.sort(corner_indices);
		for (i = 0; i < 6; i++) {
			x_corners_new[i] = x_corners[corner_indices[i]];
			y_corners_new[i] = y_corners[corner_indices[i]];			
		}
		corners_new = Array.concat(x_corners_new, y_corners_new);
	}
	return corners_new;
}

/*  if less than 6 corners are specified, find a corner between the points forming the largest angle */
// initial_corners --- array of corner coordinates [[x_coords] + [y_coords]]
function addCorners(initial_corners) {

	/*  to avoid infinite looping restrict restrict loop number to MAX_CORNER_ADD_CYCLE */
	cycle_counter = 0;			
	corners = initial_corners;
	
	while (corners.length < 12 && cycle_counter < MAX_CORNER_ADD_CYCLE) {

		x_corners = Array.slice(corners, 0, corners.length / 2);
		y_corners = Array.slice(corners, corners.length / 2);
		print(i + " adding corner found only " +  x_corners.length);
		
		/*  identify thw two points forming the largest angle with the center */
		angles = newArray(x_corners.length);
		for (c = 0; c < x_corners.length; c++) {
			angle = calculateAngle(x_corners[c], y_corners[c], x_corners[(c + 1) % x_corners.length], y_corners[(c + 1) % x_corners.length], x_c, y_c);
			angle += PI;
			angle *= 180 / PI;
			angles[c] = angle;			
		}
		Array.getStatistics(angles, min, max_angle, mean, stdDev);
		max_index = findIndex(angles, max_angle);

		/* calculate start and stop angles and try to find a new max */
		angle0 = atan2(y_corners[max_index] - y_c, x_corners[max_index] - x_c);
		angle1 = atan2(y_corners[(max_index + 1) % x_corners.length] - y_c, x_corners[(max_index + 1) % x_corners.length] - x_c);
		angle0 += PI;
		angle0 *= 180 / PI;
		angle1 += PI;
		angle1 *= 180 / PI;
		new_anglesweep = getProfilesAsArray("Particles", x_c, y_c, radius, angle0, angle1, false);
		new_max = detectMaxima(new_anglesweep, true);
						
		/* detect new edge point at that angle and insert it into the initial corner array */
		x_initial_new = x_c - radius * cos((angle0 + new_max * SAMPLING_WINDOW) * PI / 180);
		y_initial_new = y_c - radius * sin((angle0 + new_max * SAMPLING_WINDOW) * PI / 180);
		corner_new = refineCornerDistance(x_c, y_c, x_initial_new, y_initial_new);
		x_corners_new = corner_new[0];
		y_corners_new = corner_new[1];		

		x_corners_low = Array.slice(x_corners, 0, max_index + 1);
		x_corners_high = Array.slice(x_corners, max_index + 1);
		x_corners_low = Array.concat(x_corners_low, x_corners_new);
		x_corners = Array.concat(x_corners_low, x_corners_high);
		y_corners_low = Array.slice(y_corners, 0, max_index + 1);
		y_corners_high = Array.slice(y_corners, max_index + 1);
		y_corners_low = Array.concat(y_corners_low, y_corners_new);
		y_corners = Array.concat(y_corners_low, y_corners_high);
		
		corners = Array.concat(x_corners, y_corners);
		
		/* check if new point is on an edge */
		corners = removeCornersOnEdges(corners);
		
		cycle_counter++;				
	}

	return corners;	
}

/* calculate are of the hexagon from the area of the six triangular segments */
// x_c, y_c --- center coordinates
// x_corners, y_corners --- arrays containing the corner coordinates
function calculateArea(x_c, y_c, x_corners, y_corners) {

	total_area = 0;
	/* a, b, c --- side lengths of the triangle */
	for (corner = 0; corner < 6; corner++) {
		x0 = x_corners[corner];
		y0 = y_corners[corner];
		x1 = x_corners[(corner + 1) % 6];
		y1 = y_corners[(corner + 1) % 6];
		a = calculateDistance(x0, y0, x1, y1);
		b = calculateDistance(x0, y0, x_c, y_c);
		c = calculateDistance(x_c, y_c, x1, y1);
		
		/* Heron's formula */
		s = (a + b + c) / 2;		
		total_area += sqrt(s * (s - a) * (s - b) * (s - c));
	}
	
	return total_area * pow(pixel_size, 2);
}

/* calculate how much the hexagon is squeezed */
// angles --- angles between adjacent points, center as vertex
function calculateExcentricity(angles) {	
	opp_angles = newArray(3);
	for (a = 0; a < 3; a++) {
		opp_angles[a] = angles[a] + angles[(a + 3) % 6];		
	}
	Array.getStatistics(opp_angles, min_angle, max_angle, VOID, VOID);
	return min_angle / max_angle; 
}

/* calculate how much the hexagon is squeezed */
// distances --- distances between adjacent points
function calculateExcentricityfromDistances(distances) {	
	excentricities = newArray(3);
	for (r = 0; r < 3; r++) {
		two_adjacent_points = (distances[r] + distances[r + 1] + distances[(r + 3)] + distances[(r + 4) % 6]) / 4;
		next_point = (distances[r + 2] + distances[(r + 5) % 6]) / 2;
		excentricities[r] = two_adjacent_points / next_point;		
	}
	Array.getStatistics(excentricities, min_excentricity, VOID, VOID, VOID);
	return min_excentricity; 
}

/* calculate how different alternate angles are */
// angles --- angles between adjacent points, center as vertex
function calculateIrregularity(angles) {	
	angle_sum0 = angles[0] + angles[2] + angles[4];
	angle_sum1 = angles[1] + angles[3] + angles[5];
	if (angle_sum0 < angle_sum1) {
		return angle_sum0 / angle_sum1;
	}
	else {
		return angle_sum1 / angle_sum0;
	}
}

/* calculate how different alternate angles are */
// distances --- distances between adjacent points
function calculateIrregularityfromDistances(distances) {
	
	distance_sum0 = distances[0] + distances[2] + distances[4];
	distance_sum1 = distances[1] + distances[3] + distances[5];
	if (distance_sum0 < distance_sum1) {
		return distance_sum0 / distance_sum1;
	}
	else {
		return distance_sum1 / distance_sum0;
	}
}


/* Main starts here */
// GLOBAL STATICS used by functions, see above
SAMPLING_WINDOW = 5;
MAX_HITDISTANCE = 4;
DIFFERENTIATION_WINDOW = 5;
IGNORE_LIMIT = 4;
MAX_ALLOWED_ANGLE = 150;
MIN_OPP_ANGLE = 20;
MAX_CORNER_ADD_CYCLE = 5;

/* Preparation and edge detection (Laplacian) */
run("Set Measurements...", "area mean min centroid feret's redirect=None decimal=3");
image_title = getTitle();
image_path = getInfo("image.directory");
getDimensions(image_width, image_height, channels, slices, frames);
getPixelSize(pixelsize_unit, pixel_size, VOID, VOID);
run("Duplicate...", "title=LabeledImage");
run("Duplicate...", "title=Working");
run("Bandpass Filter...", "filter_large=50 filter_small=5 suppress=None tolerance=5 autoscale saturate");
run("FeatureJ Laplacian", "compute smoothing=5 detect");
setAutoThreshold("Default dark no-reset");
run("Analyze Particles...", "size=20-Infinity circularity=0.60-1.00 display clear include add");
close("Working Laplacian zero-crossings");

selectImage("LabeledImage");
run("RGB Color");

newImage("Particles", "8-bit black", image_width, image_height, 1);
setColor(255, 255, 255);
roiManager("Deselect");
roiManager("fill");

/* Save centroid and feret points for each detected hexagon */
x_centroids = newArray();
y_centroids = newArray();
x_ferets = newArray();
y_ferets = newArray();
number_particles = nResults;
for (i = 0 ; i < number_particles; i++) {
	x_centroids = Array.concat(x_centroids, getResult("X", i) / pixel_size);
	y_centroids = Array.concat(y_centroids, getResult("Y", i) / pixel_size);
	
	/* In ImageJ version 1.52n "Analyze Particles ..." returns FeretX and FeretY in pixel coordinates */
	//x_ferets = Array.concat(x_ferets, getResult("FeretX", i) / pixel_size);
	//y_ferets = Array.concat(y_ferets, getResult("FeretY", i) / pixel_size);
	x_ferets = Array.concat(x_ferets, getResult("FeretX", i));
	y_ferets = Array.concat(y_ferets, getResult("FeretY", i));		
}



/* identifying corners for each hexagon */
run("Clear Results");
result_counter = 0;

for (i = 0; i < number_particles; i++) {	
	x_f = x_ferets[i];
	y_f = y_ferets[i];
	x_c = x_centroids[i];
	y_c = y_centroids[i];
	radius = calculateDistance(x_c, y_c, x_f, y_f);
	
	/* measure the edge distance from the center for all angles */
	edge_profiles = getEdgesAsArray("Working", x_c, y_c, radius, 0, 360, true);
	
	/* measure in the masked image the intensity for all angles */
	particle_profiles = getProfilesAsArray("Particles", x_c, y_c, radius, 0, 360, true);	

	
	weighted_profiles = newArray(particle_profiles.length);
	for (v = 0; v < particle_profiles.length; v++) {
		weighted_profiles[v] = particle_profiles[v] + edge_profiles[v];
	}
			
	print("Data for image " + i);

	/* get maxima of weighted profiles --> corners */
	angle_values = detectMaxima(weighted_profiles, false);
	Array.sort(angle_values);	

	/* add a corner on the opposite side if there is none */
	angle_values = probeOppositeAngles(angle_values);

	/* calculate corner coordinates and refine them */
	corners = newArray(angle_values.length * 2);
	for (c = 0; c < angle_values.length; c++) {
		x_cor = x_c + radius * cos(angle_values[c] * SAMPLING_WINDOW * PI / 180);
		y_cor = y_c + radius * sin(angle_values[c] * SAMPLING_WINDOW * PI / 180);
		corner = refineCornerDistance(x_c, y_c, x_cor, y_cor);
		corners[c] = corner[0];
		corners[c + angle_values.length] = corner[1];		
	}

	/* create Convex Hull, remove corners inside hexagons */
	x_corners = Array.slice(corners, 0, corners.length / 2);
	y_corners = Array.slice(corners, corners.length / 2);
	makeSelection("polygon", x_corners, y_corners);
	run("Convex Hull");
	getSelectionCoordinates(x_corners, y_corners);
	corners = Array.concat(x_corners, y_corners);

	/* remove corners from edges */
	corners = removeCornersOnEdges(corners);
		
	/* if less than 6 corners are identified try to add one in the largest gap */
	if (corners.length < 12) {
		corners = addCorners(corners);					
	}

	/* annotate current hexagon */
	selectImage("LabeledImage");	
	setColor(255, 0, 0);
	fillRect(x_c - 1, y_c - 1, 3, 3);
	drawString(i, x_c + 5, y_c + 5);
		
	if (corners.length == 12) {
		setColor(0, 255, 0);	
	}
	else {
		setColor(255, 0, 0);	
	}	
	for (c = 0; c < corners.length / 2; c++) {
		fillRect(corners[c] - 1, corners[c + corners.length / 2] - 1, 3, 3);		
	}	

	/* if six corners have been found, calculate hexagon area, excentricity and irregularity */
	if (corners.length == 12) {
		point_distances = newArray(6);
		x_corners = Array.slice(corners, 0, corners.length / 2);
		y_corners = Array.slice(corners, corners.length / 2);
		for (p = 0; p < 6; p++) {
			point_distances[p] = calculateDistance(x_corners[p], y_corners[p], x_corners[(p + 1) % 6], y_corners[(p + 1) % 6]);
		}
		angles = calculateAngles(x_corners, y_corners, x_c, y_c);
		
		setResult("Hexagon #", result_counter, i);

		setResult("Centroid X", result_counter, x_c * pixel_size);
		setResult("Centroid Y", result_counter, y_c * pixel_size);

		area = calculateArea(x_c, y_c, x_corners, y_corners);
		setResult("Area", result_counter, area);
		print("Area in " + pixelsize_unit + "^2: " + area);
				
		excentricity_angles = calculateExcentricity(angles);
		setResult("Excentricity_angles", result_counter, excentricity_angles);
		excentricity_distances = calculateExcentricityfromDistances(point_distances);
		setResult("Excentricity_distances", result_counter, excentricity_distances);
		print("Hexagon excentricity (angles): " + excentricity_angles + " Hexagon excentricity (distances): " + excentricity_distances);
		
		irregularity_angles = calculateIrregularity(angles);
		setResult("Irregularity_angles", result_counter, irregularity_angles);
		irregularity_distances = calculateIrregularityfromDistances(point_distances);
		setResult("Irregularity_distances", result_counter, irregularity_distances);
		print("Hexagon irregularity (angles): " + irregularity_angles + " Hexagon irregularity (angles): " + irregularity_distances);

		Array.getStatistics(angle_values, min_angle);		
		setResult("First Angle", result_counter, min_angle * 2 * PI);

		for (v = 0; v < 6; v++) {
			setResult("Corner " + v + " X", result_counter, x_corners[v] * pixel_size);
			setResult("Corner " + v + " Y", result_counter, y_corners[v] * pixel_size);
		}
				
		print("------------------");
		result_counter++;
		updateResults();
	}		
}

/* save results as .csv and the annotated image */
saveAs("Results", image_path + File.separator  + image_title + "_Results.csv");
selectImage("LabeledImage");	
saveAs("tiff", image_path + File.separator  + image_title + "_labeled.tif");
close("Working");
close("Particles"); 