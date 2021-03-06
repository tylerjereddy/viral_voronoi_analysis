'''Library of utility functions for analyzing virus simulations with spherical Voronoi diagrams.'''
import glob
import math
import numpy
import collections
import scipy
import scipy.spatial
from scipy.spatial import SphericalVoronoi
import pandas

def produce_list_trajectories(path_string,glob_search_string):
    '''Function to produce a list of trajectories given a path string and a glob search string for the file names at that path. Include the final slash in the path_string'''
    list_of_trajectories = glob.glob(path_string + glob_search_string)
    return list_of_trajectories

def generate_random_array_spherical_generators(num_generators,sphere_radius,prng_object):
    '''Recoded using standard uniform selector over theta and acos phi, http://mathworld.wolfram.com/SpherePointPicking.html
    Same as in iPython notebook version'''
    u = prng_object.uniform(low=0,high=1,size=num_generators)
    v = prng_object.uniform(low=0,high=1,size=num_generators)
    theta_array = 2 * math.pi * u
    phi_array = numpy.arccos((2*v - 1.0))
    r_array = sphere_radius * numpy.ones((num_generators,))
    spherical_polar_data = numpy.column_stack((r_array,theta_array,phi_array))
    cartesian_random_points = convert_spherical_array_to_cartesian_array(spherical_polar_data)
    #filter out any duplicate generators:
    df_random_points = pandas.DataFrame(cartesian_random_points)
    df_random_points_no_duplicates = df_random_points.drop_duplicates()
    array_random_spherical_generators = df_random_points_no_duplicates.as_matrix()
    return array_random_spherical_generators

def convert_cartesian_array_to_spherical_array(coord_array,angle_measure='radians'):
    '''Take shape (N,3) cartesian coord_array and return an array of the same shape in spherical polar form (r, theta, phi). Based on StackOverflow response: http://stackoverflow.com/a/4116899
    use radians for the angles by default, degrees if angle_measure == 'degrees' '''
    spherical_coord_array = numpy.zeros(coord_array.shape)
    xy = coord_array[...,0]**2 + coord_array[...,1]**2
    spherical_coord_array[...,0] = numpy.sqrt(xy + coord_array[...,2]**2)
    spherical_coord_array[...,1] = numpy.arctan2(coord_array[...,1], coord_array[...,0])
    spherical_coord_array[...,2] = numpy.arccos(coord_array[...,2] / spherical_coord_array[...,0])
    if angle_measure == 'degrees':
        spherical_coord_array[...,1] = numpy.degrees(spherical_coord_array[...,1])
        spherical_coord_array[...,2] = numpy.degrees(spherical_coord_array[...,2])
    return spherical_coord_array

def convert_spherical_array_to_cartesian_array(spherical_coord_array,angle_measure='radians'):
    '''Take shape (N,3) spherical_coord_array (r,theta,phi) and return an array of the same shape in cartesian coordinate form (x,y,z). Based on the equations provided at: http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#From_spherical_coordinates
    use radians for the angles by default, degrees if angle_measure == 'degrees' '''
    cartesian_coord_array = numpy.zeros(spherical_coord_array.shape)
    #convert to radians if degrees are used in input (prior to Cartesian conversion process)
    if angle_measure == 'degrees':
        spherical_coord_array[...,1] = numpy.deg2rad(spherical_coord_array[...,1])
        spherical_coord_array[...,2] = numpy.deg2rad(spherical_coord_array[...,2])
    #now the conversion to Cartesian coords
    cartesian_coord_array[...,0] = spherical_coord_array[...,0] * numpy.cos(spherical_coord_array[...,1]) * numpy.sin(spherical_coord_array[...,2])
    cartesian_coord_array[...,1] = spherical_coord_array[...,0] * numpy.sin(spherical_coord_array[...,1]) * numpy.sin(spherical_coord_array[...,2])
    cartesian_coord_array[...,2] = spherical_coord_array[...,0] * numpy.cos(spherical_coord_array[...,2])
    return cartesian_coord_array

def calculate_haversine_distance_between_spherical_points(cartesian_array_1,cartesian_array_2,sphere_radius):
    '''Calculate the haversine-based distance between two points on the surface of a sphere. Should be more accurate than the arc cosine strategy. See, for example: http://en.wikipedia.org/wiki/Haversine_formula'''
    spherical_array_1 = convert_cartesian_array_to_spherical_array(cartesian_array_1)
    spherical_array_2 = convert_cartesian_array_to_spherical_array(cartesian_array_2)
    lambda_1 = spherical_array_1[1]
    lambda_2 = spherical_array_2[1]
    phi_1 = spherical_array_1[2]
    phi_2 = spherical_array_2[2]
    #we rewrite the standard Haversine slightly as long/lat is not the same as spherical coordinates - phi differs by pi/4
    spherical_distance = 2.0 * sphere_radius * math.asin(math.sqrt( ((1 - math.cos(phi_2-phi_1))/2.) + math.sin(phi_1) * math.sin(phi_2) * ( (1 - math.cos(lambda_2-lambda_1))/2.)  ))
    return spherical_distance

def calculate_surface_area_of_a_spherical_Voronoi_polygon(array_ordered_Voronoi_polygon_vertices,sphere_radius):
    '''Calculate the surface area of a polygon on the surface of a sphere. Based on equation provided here: http://mathworld.wolfram.com/LHuiliersTheorem.html
    Decompose into triangles, calculate excess for each'''
    #have to convert to unit sphere before applying the formula
    spherical_coordinates = convert_cartesian_array_to_spherical_array(array_ordered_Voronoi_polygon_vertices)
    spherical_coordinates[...,0] = 1.0
    array_ordered_Voronoi_polygon_vertices = convert_spherical_array_to_cartesian_array(spherical_coordinates)
    #handle nearly-degenerate vertices on the unit sphere by returning an area close to 0 -- may be better options, but this is my current solution to prevent crashes, etc.
    #seems to be relatively rare in my own work, but sufficiently common to cause crashes when iterating over large amounts of messy data
    if scipy.spatial.distance.pdist(array_ordered_Voronoi_polygon_vertices).min() < (10 ** -7):
        print 'Problematic spherical polygon SA calculation.'
        return 10 ** -8
    else:
        n = array_ordered_Voronoi_polygon_vertices.shape[0]
        #point we start from
        root_point = array_ordered_Voronoi_polygon_vertices[0]
        totalexcess = 0
        #loop from 1 to n-2, with point 2 to n-1 as other vertex of triangle
        # this could definitely be written more nicely
        b_point = array_ordered_Voronoi_polygon_vertices[1]
        root_b_dist = calculate_haversine_distance_between_spherical_points(root_point, b_point, 1.0)
        for i in 1 + numpy.arange(n - 2):
            a_point = b_point
            b_point = array_ordered_Voronoi_polygon_vertices[i+1]
            root_a_dist = root_b_dist
            root_b_dist = calculate_haversine_distance_between_spherical_points(root_point, b_point, 1.0)
            a_b_dist = calculate_haversine_distance_between_spherical_points(a_point, b_point, 1.0)
            s = (root_a_dist + root_b_dist + a_b_dist) / 2.
            totalexcess += 4 * math.atan(math.sqrt( math.tan(0.5 * s) * math.tan(0.5 * (s-root_a_dist)) * math.tan(0.5 * (s-root_b_dist)) * math.tan(0.5 * (s-a_b_dist))))
        return totalexcess * (sphere_radius ** 2)

        



class voronoi_neighbour_analysis(object):
    '''Accepts the dictionary of Voronoi diagram data structure I've been using in the ipynb and allows for parsing of neighbour properties in Voronoi diagrams.'''

    def __init__(self,voronoi_data_dict,inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays'):
        self.voronoi_data_dict = voronoi_data_dict
        self.outer_leaflet_data_key = outer_leaflet_vertex_list_key
        self.inner_leaflet_data_key = inner_leaflet_vertex_list_key

    def voronoi_cell_row_accounting(self, simplified_dict, frame_index):
        '''Produce convenience dictionary that contains the number of rows in each voronoi cell array, as this is quite costly to calculate repetitively if not abstacted. Basically, tricky to work with because the data is not homogenously structured (sub-arrays have different numbers of rows so numpy convenience functions mostly don't work).'''
        voronoi_cell_row_accounting_dict = {}
        for molecular_species_name, arrays_voronoi_cells_in_all_parsed_frames in simplified_dict.iteritems():
            array_voronoi_cell_coords_current_species_current_frame = arrays_voronoi_cells_in_all_parsed_frames[frame_index]
            current_species_Voronoi_cell_row_sizes = [voronoi_coords.shape[0] for voronoi_coords in array_voronoi_cell_coords_current_species_current_frame]
            voronoi_cell_row_accounting_dict[molecular_species_name] = current_species_Voronoi_cell_row_sizes
        return voronoi_cell_row_accounting_dict

    def identify_voronoi_neighbours(self,frame_index):
        '''Produce a data structure that associates each voronoi cell for each molecular species in each leaflet with the number of neighbours that it has in the Voronoi diagram.'''

        simplified_inner_leaflet_voronoi_data_dict = self.condense_voronoi_cell_data_by_leaflet(self.voronoi_data_dict,self.inner_leaflet_data_key)
        simplified_outer_leaflet_voronoi_data_dict = self.condense_voronoi_cell_data_by_leaflet(self.voronoi_data_dict,self.outer_leaflet_data_key)

        inner_leaflet_voronoi_row_accounting_dict = self.voronoi_cell_row_accounting(simplified_inner_leaflet_voronoi_data_dict, frame_index)            
        outer_leaflet_voronoi_row_accounting_dict = self.voronoi_cell_row_accounting(simplified_outer_leaflet_voronoi_data_dict, frame_index)            

        outer_leaflet_neighbour_dict = {}
        inner_leaflet_neighbour_dict = {}
      
        inner_leaflet_neighbour_dict = self.per_leaflet_accumulation_neighbour_data(self.inner_leaflet_data_key, self.voronoi_data_dict, inner_leaflet_neighbour_dict, simplified_inner_leaflet_voronoi_data_dict,frame_index, inner_leaflet_voronoi_row_accounting_dict)
        outer_leaflet_neighbour_dict = self.per_leaflet_accumulation_neighbour_data(self.outer_leaflet_data_key, self.voronoi_data_dict, outer_leaflet_neighbour_dict, simplified_outer_leaflet_voronoi_data_dict,frame_index, outer_leaflet_voronoi_row_accounting_dict)

        #idea is to have dictionary data structures that look like this (prototype at the moment): {'POPC': {5 : {'frequency': 20, 'list_surface_areas':[22,11,17], ...}}} where 5 is the number of neighbours
        #I think this won't handle multiple frames correctly yet, but just trying to get things started
        return (inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict)

    def per_leaflet_accumulation_neighbour_data(self, leaflet_data_key, data_dict, results_dictionary, simplified_data_dict, frame_index, voronoi_cell_row_accounting_dict):
        '''Populate results_dictionary with voronoi cell neighbour data structure for a given leaflet (and specific frame).'''

        def default_factory():
            '''For defaultdict accumulation of voronoi cell neighbour / surface area data structures.'''
            return {'frequency': 0, 'list_surface_areas': []}

        for molecular_species_name_string, subdictionary in data_dict.iteritems():
            print molecular_species_name_string, '(', leaflet_data_key, ')'
            neighbour_freq_dict = collections.defaultdict(default_factory)
            leaflet_voronoi_data_list_current_species = subdictionary[leaflet_data_key]
            list_voronoi_cells_current_frame = leaflet_voronoi_data_list_current_species[frame_index]
            #pull out the radius, which is the same for all particles because they have been projected to a perfect sphere in voronoi_analysis_loop
            radius = convert_cartesian_array_to_spherical_array(list_voronoi_cells_current_frame[0])[0,0]
            for voronoi_cell_coord_array in list_voronoi_cells_current_frame: #I'll want to find common vertices by checking all cells in current leaflet
                shape_voronoi_cell_coord_array = voronoi_cell_coord_array.shape
                dimensions_voronoi_cell_coord_array = voronoi_cell_coord_array.ndim 
                assert dimensions_voronoi_cell_coord_array == 2, "Each voronoi cell coordinate array should have two dimensions, but got {ndim}.".format(ndim = dimensions_voronoi_cell_coord_array)
                assert shape_voronoi_cell_coord_array[1] == 3, "Voronoi cell coordinates should have 3 data columns, but got {columns} data columns.".format(columns = shape_voronoi_cell_coord_array[1])

                neighbour_count_current_voronoi_cell = self.count_neighbours_current_frame(voronoi_cell_coord_array,simplified_data_dict,frame_index, voronoi_cell_row_accounting_dict)
                neighbour_freq_dict[neighbour_count_current_voronoi_cell]['frequency'] += 1 #debug: for the ZERO neighbour count, this incrementation is happening *thousands of times!!*
                surface_area_current_voronoi_cell = calculate_surface_area_of_a_spherical_Voronoi_polygon(voronoi_cell_coord_array,radius)
                neighbour_freq_dict[neighbour_count_current_voronoi_cell]['list_surface_areas'].append(surface_area_current_voronoi_cell)
            results_dictionary[molecular_species_name_string] = neighbour_freq_dict
            total_frequency = 0
            for neighbour_count, subdict in neighbour_freq_dict.iteritems():
                total_frequency += subdict['frequency']
            assert len(simplified_data_dict[molecular_species_name_string][frame_index]) == total_frequency, "Frequency mismatch for {species}.".format(species=molecular_species_name_string)
            print molecular_species_name_string, 'results dict produced.'
        print '(', leaflet_data_key, ')', 'overall results dict produced.'
        return results_dictionary

    def condense_voronoi_cell_data_by_leaflet(self, data_dict, leaflet_data_key):
        '''Produce a dictionary that contains a simple data structure for a given leaflet: {'molecular_species_name' : arrays_voronoi_cells_in_all_parsed_frames, ...}'''
        simplified_dict = {}
        for molecular_species_name_string, subdictionary in data_dict.iteritems():
            current_species_voronoi_data_list = subdictionary[leaflet_data_key]
            arrays_voronoi_cells_in_all_parsed_frames = numpy.array(current_species_voronoi_data_list) #debug: problem exists here #tricky data structure with nesting into different frame values
            simplified_dict[molecular_species_name_string] = arrays_voronoi_cells_in_all_parsed_frames
        return simplified_dict

    def count_neighbours_current_frame(self, single_voronoi_cell_array_of_coordinates, simplified_dict_for_leaflet, frame_index, voronoi_cell_row_accounting_dict):
        '''Assess the number of neighbouring Voronoi cells for the given Voronoi cell being probed in the leaflet and frame of interest.'''
        neighbour_count_current_voronoi_cell = 0
        list_arrays_all_Voronoi_cells_current_frame = [] #for raw neighbour analysis don't care about species identities, so just merge all the Voronoi cells in current frame to a single list
        list_all_voronoi_cell_row_sizes = []

        for molecular_species_name, list_arrays_voronoi_cells_in_all_parsed_frames in simplified_dict_for_leaflet.iteritems():
            list_arrays_voronoi_cell_coords_current_species_current_frame = list_arrays_voronoi_cells_in_all_parsed_frames[frame_index]
            list_arrays_all_Voronoi_cells_current_frame.extend(list_arrays_voronoi_cell_coords_current_species_current_frame)
            current_species_Voronoi_cell_row_sizes = voronoi_cell_row_accounting_dict[molecular_species_name]
            list_all_voronoi_cell_row_sizes.extend(current_species_Voronoi_cell_row_sizes)

        list_index_ranges = numpy.cumsum([0] + list_all_voronoi_cell_row_sizes) #overlapping: i.e., [0, 6, 11, 15, 20]
        list_index_range_tuples = zip(list_index_ranges[:-1],list_index_ranges[1:]) #should be i.e., [(0, 6), (6, 11), ...]
        flattened_array_all_voronoi_cell_coords_current_frame = numpy.concatenate(list_arrays_all_Voronoi_cells_current_frame)
        view_structured_array_flattened_array_voronoi_cell_coords_current_frame = flattened_array_all_voronoi_cell_coords_current_frame.view(dtype = 'f8,f8,f8').reshape(flattened_array_all_voronoi_cell_coords_current_frame.shape[0])
        single_voronoi_cell_array_of_coordinates_view = single_voronoi_cell_array_of_coordinates.view(dtype = 'f8,f8,f8').reshape(single_voronoi_cell_array_of_coordinates.shape[0])
        mask = numpy.in1d(view_structured_array_flattened_array_voronoi_cell_coords_current_frame, single_voronoi_cell_array_of_coordinates_view)
        non_zero_count_array = numpy.array([numpy.count_nonzero(mask[start:end]) for start, end in list_index_range_tuples])
        matching_vertices_current_cell = numpy.count_nonzero((non_zero_count_array > 0) & (non_zero_count_array < single_voronoi_cell_array_of_coordinates.shape[0])) #filter out exact matches to self
        neighbour_count_current_voronoi_cell = matching_vertices_current_cell

        return neighbour_count_current_voronoi_cell

    def sample_neighbour_data(self, frame_index, molecular_species_name, num_surrounding_neighbours, leaflet):
        '''Given the name of the molecular species of interest and the raw number of surrounding Voronoi neighbours, pull out the ordered polygon coordinates for the Voronoi cell of interest as well as for the neighbouring Voronoi cells. The idea is to eventually use the returned data to produce a sample plot of a region with N neighbours.''' 

        subdict = self.voronoi_data_dict[molecular_species_name]

        if leaflet == 'inner':
            key = self.inner_leaflet_data_key
        elif leaflet == 'outer':
            key = self.outer_leaflet_data_key

        list_vertex_arrays_current_molecular_species_current_leaflet_current_frame = subdict[key][frame_index]
        simplified_leaflet_voronoi_data_dict = self.condense_voronoi_cell_data_by_leaflet(self.voronoi_data_dict,key)
        leaflet_voronoi_row_accounting_dict = self.voronoi_cell_row_accounting(simplified_leaflet_voronoi_data_dict, frame_index)            

        for species_of_interest_Voronoi_cell_coords in list_vertex_arrays_current_molecular_species_current_leaflet_current_frame:
            neighbour_count_current_cell = self.count_neighbours_current_frame(species_of_interest_Voronoi_cell_coords, simplified_leaflet_voronoi_data_dict, frame_index, leaflet_voronoi_row_accounting_dict)
            #as soon as I find a Voronoi cell with the correct number of neighbours for the desired species type, just use that one
            if neighbour_count_current_cell == num_surrounding_neighbours:
                coordinates_of_interest_central_Voronoi_cell = species_of_interest_Voronoi_cell_coords
                break

        #now, given the coordinates of the central Voronoi cell of interest, identify the coordinates (and probably species names, if possible) of its neighbouring cells
        #for now, this will be a hacked version of count_neighbours_current_frame()

        def default_factory():
            return []

        dict_sample_Voronoi_cells = {'central_cell':coordinates_of_interest_central_Voronoi_cell, 'neighbours': collections.defaultdict(default_factory)} #for storing the central and neighbouring Voronoi cells

        for species_name, list_arrays_voronoi_cells_in_all_parsed_frames in simplified_leaflet_voronoi_data_dict.iteritems():
            list_arrays_voronoi_cell_coords_current_species_current_frame = list_arrays_voronoi_cells_in_all_parsed_frames[frame_index]
            current_species_Voronoi_cell_row_sizes = leaflet_voronoi_row_accounting_dict[species_name]
            list_index_ranges = numpy.cumsum([0] + current_species_Voronoi_cell_row_sizes) #overlapping: i.e., [0, 6, 11, 15, 20]
            list_index_range_tuples = zip(list_index_ranges[:-1],list_index_ranges[1:]) #should be i.e., [(0, 6), (6, 11), ...]
            flattened_array_all_voronoi_cell_coords_current_frame = numpy.concatenate(list_arrays_voronoi_cell_coords_current_species_current_frame)
            view_structured_array_flattened_array_voronoi_cell_coords_current_frame = flattened_array_all_voronoi_cell_coords_current_frame.view(dtype = 'f8,f8,f8').reshape(flattened_array_all_voronoi_cell_coords_current_frame.shape[0])
            single_voronoi_cell_array_of_coordinates_view = coordinates_of_interest_central_Voronoi_cell.view(dtype = 'f8,f8,f8').reshape(coordinates_of_interest_central_Voronoi_cell.shape[0])
            mask = numpy.in1d(view_structured_array_flattened_array_voronoi_cell_coords_current_frame, single_voronoi_cell_array_of_coordinates_view)
            non_zero_count_array = numpy.array([numpy.count_nonzero(mask[start:end]) for start, end in list_index_range_tuples])
            #I want to use non_zero_count_array to generate another mask -- the new mask should have 0 for 0 counts & counts == coordinates_of_interest_central_Voronoi_cell.shape[0] and 1 for all other counts
            mask_zeros = (non_zero_count_array == 0) & (non_zero_count_array == coordinates_of_interest_central_Voronoi_cell.shape[0])
            mask_ones = (non_zero_count_array > 0) & (non_zero_count_array < coordinates_of_interest_central_Voronoi_cell.shape[0])
            non_zero_count_array[mask_zeros] = 0
            non_zero_count_array[mask_ones] = 1
            for voronoi_cell_coord, mask_value in zip(list_arrays_voronoi_cell_coords_current_species_current_frame, non_zero_count_array):
                if mask_value == 1:
                    dict_sample_Voronoi_cells['neighbours'][species_name].append(voronoi_cell_coord)

        return dict_sample_Voronoi_cells


                




        
class voronoi_neighbour_analysis_by_type(voronoi_neighbour_analysis):
    '''Should parse the number of neighbours of *each type* of molecular species surrounding each voronoi cell in the Voronoi diagrams. This is intended as an extension of the basic functionality of the parent class, which only parses and reports the raw number of neighbours around each Voronoi cell.'''

    def count_neighbours_current_frame(self, single_voronoi_cell_array_of_coordinates, simplified_dict_for_leaflet, frame_index, voronoi_cell_row_accounting_dict):
        '''Assess the number of neighbouring Voronoi cells for the given Voronoi cell being probed in the leaflet and frame of interest.'''
        dict_neighbour_data_current_voronoi_cell = {}
        for molecular_species_name in simplified_dict_for_leaflet.keys():
            single_species_dict = {molecular_species_name : simplified_dict_for_leaflet[molecular_species_name]} #convenience data structure so I can reuse code in the parent class
            neighbour_count_current_voronoi_cell = super(voronoi_neighbour_analysis_by_type, self).count_neighbours_current_frame(single_voronoi_cell_array_of_coordinates, single_species_dict, frame_index, voronoi_cell_row_accounting_dict)
            dict_neighbour_data_current_voronoi_cell[molecular_species_name] = {'num_neighbours': neighbour_count_current_voronoi_cell} 
        return dict_neighbour_data_current_voronoi_cell

    def per_leaflet_accumulation_neighbour_data(self, leaflet_data_key, data_dict, results_dictionary, simplified_data_dict, frame_index, voronoi_cell_row_accounting_dict):
        '''Populate results_dictionary with voronoi cell neighbour data structure for a given leaflet (and specific frame).'''

        results_dictionary = {}

        for molecular_species_name_string, subdictionary in data_dict.iteritems():
            print molecular_species_name_string, '(', leaflet_data_key, ')'
            leaflet_voronoi_data_list_current_species = subdictionary[leaflet_data_key]
            list_voronoi_cells_current_frame = leaflet_voronoi_data_list_current_species[frame_index]
            radius = convert_cartesian_array_to_spherical_array(list_voronoi_cells_current_frame[0])[0,0]
            #print 'len(list_voronoi_cells_current_frame):', len(list_voronoi_cells_current_frame)
            for voronoi_cell_coord_array in list_voronoi_cells_current_frame: #I'll want to find common vertices by checking all cells in current leaflet
                shape_voronoi_cell_coord_array = voronoi_cell_coord_array.shape
                dimensions_voronoi_cell_coord_array = voronoi_cell_coord_array.ndim 
                assert dimensions_voronoi_cell_coord_array == 2, "Each voronoi cell coordinate array should have two dimensions, but got {ndim}.".format(ndim = dimensions_voronoi_cell_coord_array)
                assert shape_voronoi_cell_coord_array[1] == 3, "Voronoi cell coordinates should have 3 data columns, but got {columns} data columns.".format(columns = shape_voronoi_cell_coord_array[1])
                neighbour_count_subdictionary_by_lipid_type_current_voronoi_cell = self.count_neighbours_current_frame(voronoi_cell_coord_array,simplified_data_dict,frame_index, voronoi_cell_row_accounting_dict)
                surface_area_current_voronoi_cell = calculate_surface_area_of_a_spherical_Voronoi_polygon(voronoi_cell_coord_array, radius)
                neighbour_count_subdictionary_by_lipid_type_current_voronoi_cell['surface_area_voronoi_cell'] = surface_area_current_voronoi_cell
                #so, at this stage I have a dictionary object for a single Voronoi cell of a specific molecular species type
                #the dictionary is structured like this: {'POPC': {'num_neighbours': 0}, 'PPCE': {'num_neighbours': 0}, 'DPPE': {'num_neighbours': 3}, 'CER': {'num_neighbours': 1}, 'DUPC': {'num_neighbours': 1}, 'protein': {'num_neighbours': 0}, 'DOPS': {'num_neighbours': 0}, 'PPCS': {'num_neighbours': 2}}

                #try to accumulate the individual Voronoi cell results for a given molecular species into the results_dictionary 
                for key, value in neighbour_count_subdictionary_by_lipid_type_current_voronoi_cell.iteritems():
                    if key != 'surface_area_voronoi_cell':
                        mol_species_name = key
                        subdictionary_neighbour_data = value
                        #print 'subdictionary_neighbour_data.keys():', subdictionary_neighbour_data.keys()
                        num_neighbours_of_current_mol_species = subdictionary_neighbour_data['num_neighbours']
                        voronoi_cell_surface_area = neighbour_count_subdictionary_by_lipid_type_current_voronoi_cell['surface_area_voronoi_cell']
                        if not molecular_species_name_string in results_dictionary.keys():
                            results_dictionary[molecular_species_name_string] = {mol_species_name : {num_neighbours_of_current_mol_species: [voronoi_cell_surface_area]}}
                        else: #if there's already an entry for this molecular species type, then we need to check if there's already an entry for the neighbour of the given type
                            if not mol_species_name in results_dictionary[molecular_species_name_string].keys(): #can just initialize data structure
                                results_dictionary[molecular_species_name_string][mol_species_name] = {num_neighbours_of_current_mol_species: [voronoi_cell_surface_area]}
                            else: #now, check if there's already an entry for the num_neighbours in question
                                if not num_neighbours_of_current_mol_species in results_dictionary[molecular_species_name_string][mol_species_name].keys(): #again, can just initialize
                                    results_dictionary[molecular_species_name_string][mol_species_name][num_neighbours_of_current_mol_species] = [voronoi_cell_surface_area]
                                else: #append the surface area value for the new entry at this neighbour count
                                    list_voronoi_cell_surface_areas = results_dictionary[molecular_species_name_string][mol_species_name][num_neighbours_of_current_mol_species]
                                    list_voronoi_cell_surface_areas.append(voronoi_cell_surface_area)

            #print molecular_species_name_string, 'results dict produced.'
        #print '(', leaflet_data_key, ')', 'overall results dict produced.'
        return results_dictionary










def sum_Voronoi_cell_surface_areas(start_index,end_index,dictionary_voronoi_polygon_surface_areas):
    sum_protein_Voronoi_cell_surface_areas = 0
    sum_lipid_Voronoi_cell_surface_areas = 0
    for generator_index, Voronoi_cell_surface_area in dictionary_voronoi_polygon_surface_areas.iteritems():
        if generator_index < end_index and generator_index >= start_index: #protein Voronoi cells
            sum_protein_Voronoi_cell_surface_areas += Voronoi_cell_surface_area
        else: #all others are lipid Voronoi cells
            sum_lipid_Voronoi_cell_surface_areas += Voronoi_cell_surface_area
    return (sum_protein_Voronoi_cell_surface_areas,sum_lipid_Voronoi_cell_surface_areas)


def populate_dictionary_with_spherical_Voronoi_data(dict_Voronoi_cell_surface_areas,start_index,end_index,Voronoi_cell_surface_area_list,dict_Voronoi_polygon_vertices,list_arrays_Voronoi_cells,subdictionary_object,subdictionary_key_avg_surface_area,subdictionary_key_std_surface_area,subdictionary_key_vertex_arrays):
    for generator_index, Voronoi_cell_surface_area in dict_Voronoi_cell_surface_areas.iteritems():
        if generator_index >= start_index and generator_index < end_index:
            Voronoi_cell_surface_area_list.append(Voronoi_cell_surface_area)
            array_voronoi_cell_vertex_coords = dict_Voronoi_polygon_vertices[generator_index]
            list_arrays_Voronoi_cells.append(array_voronoi_cell_vertex_coords)
    voronoi_cell_surface_area_array = numpy.array(Voronoi_cell_surface_area_list)
    average_Voronoi_cell_surface_area = numpy.average(voronoi_cell_surface_area_array)
    std_dev_Voronoi_cell_surface_area = numpy.std(voronoi_cell_surface_area_array)
    subdictionary_object[subdictionary_key_avg_surface_area].append(average_Voronoi_cell_surface_area)
    subdictionary_object[subdictionary_key_std_surface_area].append(std_dev_Voronoi_cell_surface_area)
    subdictionary_object[subdictionary_key_vertex_arrays].append(list_arrays_Voronoi_cells) #the nested storage of coordinates like this could start consuming a lot of memory when I expand code

def calculate_surface_area_sphere(radius):
    surface_area = math.pi * 4 * (radius ** 2)
    return surface_area 


def TMD_particle_selector(input_array,molecule_type):
    '''Selects the TMD coordinate elements from the input array and combines to a simplified new numpy array with TMD particle centroid coordinates only.'''
    if molecule_type == 'lipid':
        output_array = input_array #still using the COG of the entire lipid
    elif molecule_type == 'HA': #the index numbers are based on study of topology combined with Danny's DPhil thesis
        HA_1_TMD_numpy_array = numpy.reshape(numpy.average(input_array[1110:1174],axis=0),(1,3))
        HA_2_TMD_numpy_array = numpy.reshape(numpy.average(input_array[2305:2369],axis=0),(1,3))
        HA_3_TMD_numpy_array = numpy.reshape(numpy.average(input_array[3500:3564],axis=0),(1,3))
        output_array = numpy.concatenate((HA_1_TMD_numpy_array,HA_2_TMD_numpy_array,HA_3_TMD_numpy_array))
    elif molecule_type == 'NA':
        NA_1_TMD_numpy_array = numpy.reshape(numpy.average(input_array[13:70],axis=0),(1,3))
        NA_2_TMD_numpy_array = numpy.reshape(numpy.average(input_array[1037:1094],axis=0),(1,3))
        NA_3_TMD_numpy_array = numpy.reshape(numpy.average(input_array[2059:2116],axis=0),(1,3))
        NA_4_TMD_numpy_array = numpy.reshape(numpy.average(input_array[3083:3140],axis=0),(1,3))
        output_array = numpy.concatenate((NA_1_TMD_numpy_array,NA_2_TMD_numpy_array,NA_3_TMD_numpy_array,NA_4_TMD_numpy_array))
    elif molecule_type == 'M2':
        M2_1_TMD_numpy_array = numpy.reshape(numpy.average(input_array[13:64],axis=0),(1,3))
        M2_2_TMD_numpy_array = numpy.reshape(numpy.average(input_array[106:157],axis=0),(1,3))
        M2_3_TMD_numpy_array = numpy.reshape(numpy.average(input_array[199:250],axis=0),(1,3))
        M2_4_TMD_numpy_array = numpy.reshape(numpy.average(input_array[292:343],axis=0),(1,3))
        output_array = numpy.concatenate((M2_1_TMD_numpy_array,M2_2_TMD_numpy_array,M2_3_TMD_numpy_array,M2_4_TMD_numpy_array))
    return output_array

def produce_Voronoi_area_dict(list_voronoi_polygon_vertices,estimated_sphere_radius):
    dictionary_Voronoi_region_surface_areas_for_each_generator = {}
    for generator_index, Voronoi_polygon_sorted_vertex_array in enumerate(list_voronoi_polygon_vertices):
        current_Voronoi_polygon_surface_area_on_sphere = calculate_surface_area_of_a_spherical_Voronoi_polygon(Voronoi_polygon_sorted_vertex_array,estimated_sphere_radius)
        assert current_Voronoi_polygon_surface_area_on_sphere > 0, "Obtained a surface area of zero for a Voronoi region."
        dictionary_Voronoi_region_surface_areas_for_each_generator[generator_index] = current_Voronoi_polygon_surface_area_on_sphere
    return dictionary_Voronoi_region_surface_areas_for_each_generator

def voronoi_analysis_loop(universe_object,start_frame,end_frame,skip_frame_value,PPCH_PO4_threshold=285,proteins_present='no',FORS_present='no',control_condition=None,dengue_condition=None):
    '''Generalization of the large Voronoi analysis loop so that I can easily expand my ipynb analysis to include all flu simulation replicates / conditions.'''
    #selections:
    if not dengue_condition:
        PPCH_PO4_selection = universe_object.select_atoms('resname PPCH and name PO4')
        FORS_AM2_selection = universe_object.select_atoms('resname FORS and name AM2')
        PPCH_PO4_threshold = PPCH_PO4_threshold #28.5 nm cutoff for outer leaflet assignment (see above)
        CHOL_ROH_selection = universe_object.select_atoms('resname CHOL and name ROH')
        DOPX_PO4_selection = universe_object.select_atoms('resname DOPX and name PO4')
        DOPE_PO4_selection = universe_object.select_atoms('resname DOPE and name PO4')
        POPS_PO4_selection = universe_object.select_atoms('resname POPS and name PO4')
        combined_selection_DOPE_DOPX_POPS_PO4 = DOPX_PO4_selection + DOPE_PO4_selection + POPS_PO4_selection
        all_lipid_selection = universe_object.select_atoms('resname PPCH or resname CHOL or resname POPS or resname DOPX or resname DOPE or resname FORS')

        if proteins_present == 'yes':
            all_protein_selection = universe_object.select_atoms('bynum 1:344388')
            if FORS_present == 'no':
                dictionary_headgroup_data = {'PPCH':{'selection':PPCH_PO4_selection},'CHOL':{'selection':CHOL_ROH_selection},'DOPX':{'selection':DOPX_PO4_selection},'DOPE':{'selection':DOPE_PO4_selection},'POPS':{'selection':POPS_PO4_selection},'protein':{'selection':all_protein_selection}}
            else:
                dictionary_headgroup_data = {'PPCH':{'selection':PPCH_PO4_selection},'CHOL':{'selection':CHOL_ROH_selection},'DOPX':{'selection':DOPX_PO4_selection},'DOPE':{'selection':DOPE_PO4_selection},'POPS':{'selection':POPS_PO4_selection},'protein':{'selection':all_protein_selection},'FORS':{'selection':FORS_AM2_selection}}
            #for virion simulation I want to assess the amount of surface area occupied by protein and lipid separately (can always sum together later to get the overall total)
            list_percent_surface_area_reconstitution_from_lipids_only = [] #outer leaflet
            list_percent_surface_area_reconstitution_from_proteins_only = [] #outer leaflet
            list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet = []
            list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet = []
        else:
            dictionary_headgroup_data = {'PPCH':{'selection':PPCH_PO4_selection},'CHOL':{'selection':CHOL_ROH_selection},'DOPX':{'selection':DOPX_PO4_selection},'DOPE':{'selection':DOPE_PO4_selection},'POPS':{'selection':POPS_PO4_selection}}
            list_percent_surface_area_reconstitution = []
            list_percent_surface_area_reconstitution_inner_leaflet = []
    else: #selections and initial data structures for dengue case
        POPC_PO4_selection = universe_object.select_atoms('resname POPC and name PO4')
        PPCE_PO4_selection = universe_object.select_atoms('resname PPCE and name PO4')
        DPPE_PO4_selection = universe_object.select_atoms('resname DPPE and name PO4')
        CER_AM2_selection = universe_object.select_atoms('resname CER and name AM2') #may have to treat ceramide differently without PO4 in headgroup region?
        DUPC_PO4_selection = universe_object.select_atoms('resname DUPC and name PO4')
        DOPS_PO4_selection = universe_object.select_atoms('resname DOPS and name PO4')
        PPCS_PO4_selection = universe_object.select_atoms('resname PPCS and name PO4')
        combined_dengue_lipid_headgroup_selection = POPC_PO4_selection + PPCE_PO4_selection + DPPE_PO4_selection + CER_AM2_selection + DUPC_PO4_selection + DOPS_PO4_selection + PPCS_PO4_selection
        all_lipid_selection = universe_object.select_atoms('resname POPC or resname PPCE or resname DPPE or resname CER or resname DUPC or resname DOPS or resname PPCS')
        all_protein_selection = universe_object.select_atoms('bynum 1:221040') #I think I'm actually going to need to select the TMDs separately (I should have some previous code written for this task -- may be able to split into TMD selections downstream, or may be easier to do it here-we'll see)
        dictionary_headgroup_data = {'POPC':{'selection':POPC_PO4_selection},'PPCE':{'selection':PPCE_PO4_selection},'DPPE':{'selection':DPPE_PO4_selection},'CER':{'selection':CER_AM2_selection},'DUPC':{'selection':DUPC_PO4_selection},'DOPS':{'selection':DOPS_PO4_selection},'PPCS':{'selection':PPCS_PO4_selection},'protein':{'selection':all_protein_selection}}
        list_percent_surface_area_reconstitution_from_lipids_only = [] #outer leaflet
        list_percent_surface_area_reconstitution_from_proteins_only = [] #outer leaflet
        list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet = []
        list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet = []

    #set up preliminary data structures
    for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
        subdictionary['voronoi_cell_avg_values_list'] = []
        subdictionary['voronoi_cell_std_values_list'] = []
        subdictionary['voronoi_cell_list_vertex_arrays'] = []
        #and now for inner leaflet data structures:
        subdictionary['voronoi_cell_avg_values_list_inner_leaflet'] = []
        subdictionary['voronoi_cell_std_values_list_inner_leaflet'] = []
        subdictionary['voronoi_cell_list_vertex_arrays_inner_leaflet'] = []

    list_frame_numbers = []

    simulation_trajectory_object = universe_object.trajectory
    if end_frame == 'full':
        end_frame = int(simulation_trajectory_object.n_frames)
    for ts in simulation_trajectory_object[start_frame:end_frame:skip_frame_value]: 
        lipid_centroid = all_lipid_selection.centroid()
        if not dengue_condition:
            PPCH_PO4_coords = dictionary_headgroup_data['PPCH']['selection'].coordinates() - lipid_centroid
            PPCH_PO4_spherical_coords = convert_cartesian_array_to_spherical_array(PPCH_PO4_coords)
        if not control_condition and not dengue_condition:
            outer_leaflet_projection_radius = numpy.average(PPCH_PO4_spherical_coords[...,0])
        elif dengue_condition:
            outer_leaflet_projection_radius = 214.6 #I calculated this as the average max dengue radial lipid headgroup distance (across all simulation frames)
        else:
            outer_leaflet_projection_radius = 700.

        if not dengue_condition:
            #use the same inner leaflet projection criterion that was used for sim33
            combined_DOPE_DOPX_POPS_PO4_coords = combined_selection_DOPE_DOPX_POPS_PO4.coordinates() - lipid_centroid
            combined_DOPE_DOPX_POPS_PO4_spherical_coords = convert_cartesian_array_to_spherical_array(combined_DOPE_DOPX_POPS_PO4_coords)
            max_DOPE_DOPX_POPS_PO4_radial_distance = numpy.sort(combined_DOPE_DOPX_POPS_PO4_spherical_coords[...,0])[-2]
            min_DOPE_DOPX_POPS_PO4_radial_distance = combined_DOPE_DOPX_POPS_PO4_spherical_coords[...,0].min()
            unbiased_midpoint_radial_distance_DOPE_DOPX_POPS_PO4 = (max_DOPE_DOPX_POPS_PO4_radial_distance + min_DOPE_DOPX_POPS_PO4_radial_distance) / 2.
            inner_leaflet_DOPE_DOPX_POPS_PO4_spherical_coords = combined_DOPE_DOPX_POPS_PO4_spherical_coords[combined_DOPE_DOPX_POPS_PO4_spherical_coords[...,0] < unbiased_midpoint_radial_distance_DOPE_DOPX_POPS_PO4]

        if not control_condition and not dengue_condition:
            inner_leaflet_projection_radius = numpy.average(inner_leaflet_DOPE_DOPX_POPS_PO4_spherical_coords[...,0])
        elif dengue_condition:
            inner_leaflet_projection_radius = 145.7 #calculated as the average minimum lipid headgroup radial distance across all dengue simulation frames
        else:
            inner_leaflet_projection_radius = 500. 
        index_counter = 0
        inner_leaflet_index_counter = 0
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            current_headgroup_MDA_selection = subdictionary['selection'] #protein selection is hardly 'headgroup,' but same treatment
            assert current_headgroup_MDA_selection.n_atoms > 0, "Number of selected {resname} headgroup particles not greater than 0.".format(resname=residue_name)
            current_headgroup_coordinate_array = current_headgroup_MDA_selection.coordinates()
            current_headgroup_coordinate_array -= lipid_centroid #center at origin
            current_headgroup_spherical_polar_coord_array = convert_cartesian_array_to_spherical_array(current_headgroup_coordinate_array)
            #perform the necessary filtering and projection based on residue type
            if residue_name == 'PPCH' or residue_name == 'FORS' and not dengue_condition:
                if control_condition:
                    PPCH_PO4_threshold = 600.
                outer_leaflet_spherical_coord_array = current_headgroup_spherical_polar_coord_array[current_headgroup_spherical_polar_coord_array[...,0] > PPCH_PO4_threshold]
                outer_leaflet_spherical_coord_array[...,0] = outer_leaflet_projection_radius
                inner_leaflet_spherical_coord_array = current_headgroup_spherical_polar_coord_array[current_headgroup_spherical_polar_coord_array[...,0] < PPCH_PO4_threshold]
                inner_leaflet_spherical_coord_array[...,0] = inner_leaflet_projection_radius
            elif residue_name == 'protein' and not dengue_condition: #now trying a strategy that isolates the protein TMDs for projection
                HA_all_particles_coord_array = current_headgroup_coordinate_array[0:289680,...]
                NA_all_particles_coord_array = current_headgroup_coordinate_array[289680:338808,...]
                M2_all_particles_coord_array = current_headgroup_coordinate_array[338808:344388,...]
                list_individual_HA_protein_coordinate_arrays = numpy.split(HA_all_particles_coord_array,80) #split to list of coord arrays for each of the 80 HA molecules
                list_individual_NA_protein_coordinate_arrays = numpy.split(NA_all_particles_coord_array,12) #similarly for NA
                list_individual_M2_protein_coordinate_arrays = numpy.split(M2_all_particles_coord_array,15) #similarly for M2 
                list_HA_TMD_coordinate_arrays = [TMD_particle_selector(HA_coord_array,'HA') for HA_coord_array in list_individual_HA_protein_coordinate_arrays]
                list_NA_TMD_coordinate_arrays = [TMD_particle_selector(NA_coord_array,'NA') for NA_coord_array in list_individual_NA_protein_coordinate_arrays]
                list_M2_TMD_coordinate_arrays = [TMD_particle_selector(M2_coord_array,'M2') for M2_coord_array in list_individual_M2_protein_coordinate_arrays] 
                array_HA_TMD_centroids = list_HA_TMD_coordinate_arrays[0]
                for HA_TMD_array in list_HA_TMD_coordinate_arrays[1:]:
                    array_HA_TMD_centroids = numpy.concatenate((array_HA_TMD_centroids,HA_TMD_array))
                assert array_HA_TMD_centroids.shape == (240,3), "There should be 240 HA TMD controids, from 80 HA x 3 TMD per trimer, but got shape {shape} and debug array shape {debug_shape}.".format(shape=array_HA_TMD_centroids.shape, debug_shape = debug_array_HA.shape)
                array_NA_TMD_centroids = list_NA_TMD_coordinate_arrays[0]
                for NA_TMD_array in list_NA_TMD_coordinate_arrays[1:]:
                    array_NA_TMD_centroids = numpy.concatenate((array_NA_TMD_centroids,NA_TMD_array))
                assert array_NA_TMD_centroids.shape == (48,3), "There should be 48 NA TMD controids, from 12 NA x 4 TMD per tetramer., but got shape {shape}".format(shape = array_NA_TMD_centroids.shape)
                array_M2_TMD_centroids = list_M2_TMD_coordinate_arrays[0]
                for M2_TMD_array in list_M2_TMD_coordinate_arrays[1:]:
                    array_M2_TMD_centroids = numpy.concatenate((array_M2_TMD_centroids,M2_TMD_array))
                assert array_M2_TMD_centroids.shape == (60,3), "There should be 60 M2 centroids, from 15 M2 x 4 TMD per tetramer, but got shape {shape}.".format(shape = array_M2_TMD_centroids.shape)
                #concatenate HA, NA, M2 TMD centroid coords to a single array (using previous nomenclature)
                current_headgroup_coordinate_array = numpy.concatenate((array_HA_TMD_centroids,array_NA_TMD_centroids,array_M2_TMD_centroids))
                assert current_headgroup_coordinate_array.shape == (348,3), "There should be 348 centroids for HA, NA and M2 TMDs in 3 dimensions." 
                current_headgroup_spherical_polar_coord_array = convert_cartesian_array_to_spherical_array(current_headgroup_coordinate_array)
                #crudely project the TMD centroid particles both up AND down to fill in proteins spaces in both leaflets (how 'bad' is this for area approximation in Voronoi diagram?!)
                outer_leaflet_spherical_coord_array = numpy.copy(current_headgroup_spherical_polar_coord_array)
                outer_leaflet_spherical_coord_array[...,0] = outer_leaflet_projection_radius
                inner_leaflet_spherical_coord_array = numpy.copy(current_headgroup_spherical_polar_coord_array)
                inner_leaflet_spherical_coord_array[...,0] = inner_leaflet_projection_radius
            elif residue_name == 'protein' and dengue_condition: #isolate dengue protein TMD centroids most likely
                list_asymmetric_unit_coordinate_arrays = numpy.split(current_headgroup_coordinate_array,60) #split into the 60 asymmetric units
                list_asymmetric_unit_EM_TMD_centroid_lists = [TMD_particle_selector_dengue(asym_unit_coords) for asym_unit_coords in list_asymmetric_unit_coordinate_arrays]
                list_EM_TMD_centroids_combined = []
                for EM_TMD_centroid_list in list_asymmetric_unit_EM_TMD_centroid_lists:
                    list_EM_TMD_centroids_combined.extend(EM_TMD_centroid_list[0])
                    list_EM_TMD_centroids_combined.extend(EM_TMD_centroid_list[1])
                current_headgroup_coordinate_array = numpy.array(list_EM_TMD_centroids_combined)
                assert current_headgroup_coordinate_array.shape == (720,3), "The dengue virion should have 720 TMD centroids."
                #crudely project the TMD centroids to both leaflets
                current_headgroup_spherical_polar_coord_array = convert_cartesian_array_to_spherical_array(current_headgroup_coordinate_array)
                outer_leaflet_spherical_coord_array = numpy.copy(current_headgroup_spherical_polar_coord_array)
                outer_leaflet_spherical_coord_array[...,0] = outer_leaflet_projection_radius
                inner_leaflet_spherical_coord_array = numpy.copy(current_headgroup_spherical_polar_coord_array)
                inner_leaflet_spherical_coord_array[...,0] = inner_leaflet_projection_radius
            else: #all other residues use a midpoint filtering method (I think I can use as-is for dengue as well?!)
                sorted_radial_distance_array = numpy.sort(current_headgroup_spherical_polar_coord_array[...,0])
                conservative_max_value = sorted_radial_distance_array[-2] #avoid floater
                conservative_min_value = sorted_radial_distance_array[1] #avoid floater
                midpoint_radial_distance = (conservative_min_value + conservative_max_value) / 2.
                if control_condition:
                    midpoint_radial_distance = 600.
                outer_leaflet_spherical_coord_array = current_headgroup_spherical_polar_coord_array[current_headgroup_spherical_polar_coord_array[...,0] > midpoint_radial_distance]
                outer_leaflet_spherical_coord_array[...,0] = outer_leaflet_projection_radius
                inner_leaflet_spherical_coord_array = current_headgroup_spherical_polar_coord_array[current_headgroup_spherical_polar_coord_array[...,0] < midpoint_radial_distance]
                inner_leaflet_spherical_coord_array[...,0] = inner_leaflet_projection_radius
        
            if index_counter == 0: #initialize the projected outer leaflet coord array rather than concatenating for first residue type
                projected_outer_leaflet_coordinate_array = convert_spherical_array_to_cartesian_array(outer_leaflet_spherical_coord_array)
            else:
                projected_outer_leaflet_coordinate_array = numpy.concatenate((projected_outer_leaflet_coordinate_array,convert_spherical_array_to_cartesian_array(outer_leaflet_spherical_coord_array)))
            #also need to track the coordinate indices for data structure management with Voronoi code (i.e., which cell areas correspond to which residue types)
            dictionary_headgroup_data[residue_name]['start_index'] = index_counter
            dictionary_headgroup_data[residue_name]['end_index'] = index_counter + outer_leaflet_spherical_coord_array.shape[0]
            index_counter += outer_leaflet_spherical_coord_array.shape[0]
            
            
            if inner_leaflet_index_counter == 0: 
                projected_inner_leaflet_coordinate_array = convert_spherical_array_to_cartesian_array(inner_leaflet_spherical_coord_array)
            else:
                projected_inner_leaflet_coordinate_array = numpy.concatenate((projected_inner_leaflet_coordinate_array,convert_spherical_array_to_cartesian_array(inner_leaflet_spherical_coord_array)))
            dictionary_headgroup_data[residue_name]['inner_leaflet_start_index'] = inner_leaflet_index_counter
            dictionary_headgroup_data[residue_name]['inner_leaflet_end_index'] = inner_leaflet_index_counter + inner_leaflet_spherical_coord_array.shape[0]
            inner_leaflet_index_counter += inner_leaflet_spherical_coord_array.shape[0]
        
        voronoi_instance = SphericalVoronoi(projected_outer_leaflet_coordinate_array,outer_leaflet_projection_radius)
        inner_leaflet_voronoi_instance = SphericalVoronoi(projected_inner_leaflet_coordinate_array,inner_leaflet_projection_radius)
        voronoi_instance.sort_vertices_of_regions()
        inner_leaflet_voronoi_instance.sort_vertices_of_regions()
        list_voronoi_polygon_vertex_indices = voronoi_instance.regions #for sample plotting Voronoi diagrams
        list_voronoi_polygon_vertex_indices_inner_leaflet = inner_leaflet_voronoi_instance.regions #for sample plotting Voronoi diagrams
        list_voronoi_polygon_vertices = []
        list_voronoi_polygon_vertices_inner_leaflet = []
        for region in list_voronoi_polygon_vertex_indices:
            list_voronoi_polygon_vertices.append(voronoi_instance.vertices[region])
        for region in list_voronoi_polygon_vertex_indices_inner_leaflet:
            list_voronoi_polygon_vertices_inner_leaflet.append(inner_leaflet_voronoi_instance.vertices[region])

        dictionary_voronoi_polygon_surface_areas = produce_Voronoi_area_dict(list_voronoi_polygon_vertices,voronoi_instance.radius)
        dictionary_voronoi_polygon_surface_areas_inner_leaflet = produce_Voronoi_area_dict(list_voronoi_polygon_vertices_inner_leaflet,inner_leaflet_voronoi_instance.radius)
        frame_number = ts.frame
        theoretical_surface_area = calculate_surface_area_sphere(outer_leaflet_projection_radius)
        theoretical_surface_area_inner_leaflet = calculate_surface_area_sphere(inner_leaflet_projection_radius)
        #for the virion I'll want to assess % SA reconstitution separately for lipid and protein
        #need to determine which generator indices correspond to protein vs lipid for this splitting
        if proteins_present == 'yes':
            protein_outer_leaflet_start_index = dictionary_headgroup_data['protein']['start_index'] 
            protein_outer_leaflet_end_index = dictionary_headgroup_data['protein']['end_index'] 
            sum_protein_Voronoi_cell_surface_areas,sum_lipid_Voronoi_cell_surface_areas = sum_Voronoi_cell_surface_areas(protein_outer_leaflet_start_index,protein_outer_leaflet_end_index,dictionary_voronoi_polygon_surface_areas)

            protein_inner_leaflet_start_index = dictionary_headgroup_data['protein']['inner_leaflet_start_index'] 
            protein_inner_leaflet_end_index = dictionary_headgroup_data['protein']['inner_leaflet_end_index'] 
            sum_protein_Voronoi_cell_surface_areas_inner_leaflet, sum_lipid_Voronoi_cell_surface_areas_inner_leaflet = sum_Voronoi_cell_surface_areas(protein_inner_leaflet_start_index,protein_inner_leaflet_end_index,dictionary_voronoi_polygon_surface_areas_inner_leaflet)

                    
            percent_surface_area_reconstitution_proteins_only = (sum_protein_Voronoi_cell_surface_areas / theoretical_surface_area) * 100
            percent_surface_area_reconstitution_proteins_only_inner_leaflet = (sum_protein_Voronoi_cell_surface_areas_inner_leaflet / theoretical_surface_area_inner_leaflet) * 100
            list_percent_surface_area_reconstitution_from_proteins_only.append(percent_surface_area_reconstitution_proteins_only)
            list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet.append(percent_surface_area_reconstitution_proteins_only_inner_leaflet)
            
            percent_surface_area_reconstitution_lipids_only = (sum_lipid_Voronoi_cell_surface_areas / theoretical_surface_area) * 100
            percent_surface_area_reconstitution_lipids_only_inner_leaflet = (sum_lipid_Voronoi_cell_surface_areas_inner_leaflet / theoretical_surface_area_inner_leaflet) * 100
            
            list_percent_surface_area_reconstitution_from_lipids_only.append(percent_surface_area_reconstitution_lipids_only)
            list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet.append(percent_surface_area_reconstitution_lipids_only_inner_leaflet)
        else:
            percent_surface_area_reconstitution = (sum(dictionary_voronoi_polygon_surface_areas.itervalues()) / theoretical_surface_area) * 100.
            percent_surface_area_reconstitution_inner_leaflet = (sum(dictionary_voronoi_polygon_surface_areas_inner_leaflet.itervalues()) / theoretical_surface_area_inner_leaflet) * 100.
            list_percent_surface_area_reconstitution.append(percent_surface_area_reconstitution)
            list_percent_surface_area_reconstitution_inner_leaflet.append(percent_surface_area_reconstitution_inner_leaflet)

        list_frame_numbers.append(frame_number)
        
        #now try to slot the Voronoi cell surface area data into the appropriate dictionary entries based on index values
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            #print 'residue_name:', residue_name
            voronoi_cell_surface_area_list = [] 
            voronoi_cell_surface_area_list_inner_leaflet = [] 
            list_arrays_Voronoi_cells = [] #for plotting purposes
            list_arrays_Voronoi_cells_inner_leaflet = [] #for plotting purposes
            start_index = subdictionary['start_index']
            end_index = subdictionary['end_index']
            inner_leaflet_start_index = subdictionary['inner_leaflet_start_index']
            inner_leaflet_end_index = subdictionary['inner_leaflet_end_index']
            #outer leaflet:
            #print 'populating outer leaflet dict data:'
            populate_dictionary_with_spherical_Voronoi_data(dict_Voronoi_cell_surface_areas=dictionary_voronoi_polygon_surface_areas,start_index=start_index,end_index=end_index,Voronoi_cell_surface_area_list=voronoi_cell_surface_area_list,dict_Voronoi_polygon_vertices=list_voronoi_polygon_vertices,list_arrays_Voronoi_cells=list_arrays_Voronoi_cells,subdictionary_object=subdictionary,subdictionary_key_avg_surface_area='voronoi_cell_avg_values_list',subdictionary_key_std_surface_area='voronoi_cell_std_values_list',subdictionary_key_vertex_arrays='voronoi_cell_list_vertex_arrays')
            #inner leaflet:
            #print 'populating inner leaflet dict data:'
            populate_dictionary_with_spherical_Voronoi_data(dict_Voronoi_cell_surface_areas=dictionary_voronoi_polygon_surface_areas_inner_leaflet,start_index=inner_leaflet_start_index,end_index=inner_leaflet_end_index,Voronoi_cell_surface_area_list=voronoi_cell_surface_area_list_inner_leaflet,dict_Voronoi_polygon_vertices=list_voronoi_polygon_vertices_inner_leaflet,list_arrays_Voronoi_cells=list_arrays_Voronoi_cells_inner_leaflet,subdictionary_object=subdictionary,subdictionary_key_avg_surface_area='voronoi_cell_avg_values_list_inner_leaflet',subdictionary_key_std_surface_area='voronoi_cell_std_values_list_inner_leaflet',subdictionary_key_vertex_arrays='voronoi_cell_list_vertex_arrays_inner_leaflet')
            
        print 'frame:', frame_number
    residue_name_list = dictionary_headgroup_data.keys()
    for key in residue_name_list:
        del dictionary_headgroup_data[key]['selection']

    if proteins_present == 'yes':
        return (list_frame_numbers,list_percent_surface_area_reconstitution_from_lipids_only,list_percent_surface_area_reconstitution_from_proteins_only,list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet,list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet,dictionary_headgroup_data)
    else:
        return (list_frame_numbers,list_percent_surface_area_reconstitution,list_percent_surface_area_reconstitution_inner_leaflet,dictionary_headgroup_data)


def produce_universe_object_on_remote_engine(data_path_1 = None,data_path_2 = None,limit_1=None,limit_2=None,limit_3=None,limit_4 = None,coordinate_filepath=None,traj_data_extension=None,traj_data_extension_with_replace=None):
    '''For loading MDA universe object on a remote engine.'''
    import MDAnalysis
    import numpy
    import scipy
    import math 
    #produce a list of trajectory files:
    list_trajectories_compact_no_solvent = sorted(produce_list_trajectories(data_path_1,'*no_solvent*xtc'),key= lambda file_string: int(file_string[limit_1:limit_2].replace('_',''))) #sort by file name part number
    if traj_data_extension: #extend the list of trajectories, if applicable
        list_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories(data_path_2,'*no_solvent*xtc'),key= lambda file_string: int(file_string[limit_3:limit_4])))
    elif traj_data_extension_with_replace: #in some cases have to use the different format with replacement method
        list_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories(data_path_2,'*no_solvent*xtc'),key= lambda file_string: int(file_string[limit_3:limit_4].replace('_',''))))
    universe_object = MDAnalysis.Universe(coordinate_filepath,list_trajectories_compact_no_solvent) 
    return universe_object

def produce_universe_object_on_remote_engine_dengue(coordinate_file_path, traj_file_path = '/sansom/n22/bioc1009/sim126_extended/translated_test.xtc'):
    '''Produce MDA Universe object on a remote IPython engine for dengue simulation [sim126 extended].'''
    import MDAnalysis
    import numpy
    import scipy
    import math 
    universe_object = MDAnalysis.Universe(coordinate_file_path, traj_file_path)
    return universe_object


def precursor_radial_distance_analysis_dengue(universe_object):
    '''Modified version of precursor_radial_distance_analysis() intended for analysis of dengue virion simulation. I think all dengue lipids should be assessed in the same manner because a symmetrical lipid species distribution (ER-derived) was assumed during the construction process.'''
    import MDAnalysis.core.distances
    POPC_PO4_selection = universe_object.select_atoms('resname POPC and name PO4')
    PPCE_PO4_selection = universe_object.select_atoms('resname PPCE and name PO4')
    DPPE_PO4_selection = universe_object.select_atoms('resname DPPE and name PO4')
    CER_AM2_selection = universe_object.select_atoms('resname CER and name AM2') #may have to treat ceramide differently without PO4 in headgroup region?
    DUPC_PO4_selection = universe_object.select_atoms('resname DUPC and name PO4')
    DOPS_PO4_selection = universe_object.select_atoms('resname DOPS and name PO4')
    PPCS_PO4_selection = universe_object.select_atoms('resname PPCS and name PO4')
    combined_dengue_lipid_selection = POPC_PO4_selection + PPCE_PO4_selection + DPPE_PO4_selection + CER_AM2_selection + DUPC_PO4_selection + DOPS_PO4_selection + PPCS_PO4_selection
    protein_selection = universe_object.select_atoms('bynum 1:221040')
    total_dengue_particles_assessed = combined_dengue_lipid_selection.n_atoms
    threshold = 275 #adjust as needed to split leaflet populations appropriately
    all_lipid_selection = universe_object.select_atoms('resname POPC or resname PPCE or resname DPPE or resname CER or resname DUPC or resname DOPS or resname PPCS')
   
    list_min_protein_distances = []
    list_max_protein_distances = []
    list_min_dengue_lipid_headgroup_distances = []
    list_max_dengue_lipid_headgroup_distances = []
    list_average_dengue_lipid_headgroup_distances = []
    list_std_dev_dengue_lipid_headgroup_distances = []
    list_dengue_lipid_headgroup_midpoint_distances = []
    list_dengue_lipid_headgroup_percent_above_threshold = []
    list_dengue_lipid_headgroup_percent_below_threshold = []
    list_frame_numbers = []
    #debug_coords_written = 0

    for ts in universe_object.trajectory[::10]: #every 10th frame (roughly every 10 ns for skip-10-filtered dengue trajectory)
        dengue_lipid_headgroup_coordinates = combined_dengue_lipid_selection.coordinates()
        dengue_protein_coordinates = protein_selection.coordinates()
        all_lipid_centroid = all_lipid_selection.centroid()
        #place the centroid of the system at the origin
        dengue_lipid_headgroup_coordinates -= all_lipid_centroid
        dengue_protein_coordinates -= all_lipid_centroid
        spherical_polar_dengue_lipid_headgroup_coordinates = convert_cartesian_array_to_spherical_array(dengue_lipid_headgroup_coordinates)
        spherical_polar_dengue_protein_coordinates = convert_cartesian_array_to_spherical_array(dengue_protein_coordinates)
        #assess the positions of the dengue lipid headgroup particles
        minimum_dengue_lipid_headgroup_radial_distance = spherical_polar_dengue_lipid_headgroup_coordinates[...,0].min()
        maximum_dengue_lipid_headgroup_radial_distance = numpy.sort(spherical_polar_dengue_lipid_headgroup_coordinates[...,0])[-1] #looks like we have a DUPC floater, based on visual inspection of debug coords printed below
        min_protein_distance = spherical_polar_dengue_protein_coordinates[...,0].min()
        max_protein_distance = spherical_polar_dengue_protein_coordinates[...,0].max()
        #debug possible floater(s) at unusually large radial distances:
        #if debug_coords_written < 1 and maximum_dengue_lipid_headgroup_radial_distance > 310:
            #import MDAnalysis.coordinates.GRO
            #writer_instance = MDAnalysis.coordinates.GRO.GROWriter('/sansom/n22/bioc1009/spherical_Voronoi_virus_work/dengue_large_radial_distance_debug.gro')
            #writer_instance.write(combined_dengue_lipid_selection)
            #debug_coords_written += 1 #only want one debug coord set for now
        midpoint_dengue_lipid_headgroup_radial_distance = numpy.average(numpy.array([minimum_dengue_lipid_headgroup_radial_distance,maximum_dengue_lipid_headgroup_radial_distance])) #unbiased midpoint for upper / lower leaflet cutoff
        average_dengue_lipid_headgroup_radial_distance = numpy.average(spherical_polar_dengue_lipid_headgroup_coordinates[...,0])
        std_dev_dengue_lipid_headgroup_radial_distance = numpy.std(spherical_polar_dengue_lipid_headgroup_coordinates[...,0])
        number_of_dengue_lipid_headgroup_radial_distances_above_threshold = (spherical_polar_dengue_lipid_headgroup_coordinates[...,0] > midpoint_dengue_lipid_headgroup_radial_distance).sum()
        number_of_dengue_lipid_headgroup_radial_distances_below_threshold = (spherical_polar_dengue_lipid_headgroup_coordinates[...,0] < midpoint_dengue_lipid_headgroup_radial_distance).sum()
        percent_dengue_lipid_headgroup_above_threshold = float(number_of_dengue_lipid_headgroup_radial_distances_above_threshold) / float(total_dengue_particles_assessed) * 100.
        percent_dengue_lipid_headgroup_below_threshold = float(number_of_dengue_lipid_headgroup_radial_distances_below_threshold) / float(total_dengue_particles_assessed) * 100.
        list_min_dengue_lipid_headgroup_distances.append(minimum_dengue_lipid_headgroup_radial_distance)
        list_max_dengue_lipid_headgroup_distances.append(maximum_dengue_lipid_headgroup_radial_distance)
        list_average_dengue_lipid_headgroup_distances.append(average_dengue_lipid_headgroup_radial_distance)
        list_std_dev_dengue_lipid_headgroup_distances.append(std_dev_dengue_lipid_headgroup_radial_distance)
        list_dengue_lipid_headgroup_midpoint_distances.append(midpoint_dengue_lipid_headgroup_radial_distance)
        list_dengue_lipid_headgroup_percent_above_threshold.append(percent_dengue_lipid_headgroup_above_threshold)
        list_dengue_lipid_headgroup_percent_below_threshold.append(percent_dengue_lipid_headgroup_below_threshold)
        list_min_protein_distances.append(min_protein_distance)
        list_max_protein_distances.append(max_protein_distance)

        frame_number = ts.frame
        list_frame_numbers.append(frame_number)

    array_min_dengue_lipid_headgroup_radial_distances = numpy.array(list_min_dengue_lipid_headgroup_distances)/ 10.
    array_max_dengue_lipid_headgroup_radial_distances = numpy.array(list_max_dengue_lipid_headgroup_distances)/ 10.
    array_min_protein_distances = numpy.array(list_min_protein_distances) / 10.
    array_max_protein_distances = numpy.array(list_max_protein_distances) / 10.
    array_average_dengue_lipid_headgroup_radial_distances = numpy.array(list_average_dengue_lipid_headgroup_distances)/ 10.
    array_std_dev_dengue_lipid_headgroup_radial_distances = numpy.array(list_std_dev_dengue_lipid_headgroup_distances)/ 10.
    array_dengue_lipid_headgroup_unbiased_midpoint_distances = numpy.array(list_dengue_lipid_headgroup_midpoint_distances) / 10.
    array_dengue_lipid_headgroup_percent_above_midpoint_threshold = numpy.array(list_dengue_lipid_headgroup_percent_above_threshold)
    array_dengue_lipid_headgroup_percent_below_midpoint_threshold = numpy.array(list_dengue_lipid_headgroup_percent_below_threshold)
    array_frame_numbers = numpy.array(list_frame_numbers)

    return (list_min_dengue_lipid_headgroup_distances,list_max_dengue_lipid_headgroup_distances,list_average_dengue_lipid_headgroup_distances,list_std_dev_dengue_lipid_headgroup_distances,list_frame_numbers,list_dengue_lipid_headgroup_percent_above_threshold,list_dengue_lipid_headgroup_percent_below_threshold,list_dengue_lipid_headgroup_midpoint_distances,list_min_protein_distances,list_max_protein_distances, array_min_dengue_lipid_headgroup_radial_distances,array_max_dengue_lipid_headgroup_radial_distances,array_min_protein_distances,array_max_protein_distances,array_average_dengue_lipid_headgroup_radial_distances,array_std_dev_dengue_lipid_headgroup_radial_distances,array_dengue_lipid_headgroup_unbiased_midpoint_distances,array_dengue_lipid_headgroup_percent_above_midpoint_threshold,array_dengue_lipid_headgroup_percent_below_midpoint_threshold,array_frame_numbers)



    



def precursor_radial_distance_analysis(universe_object,FORS_present=None,control_condition=None):
    '''This function should parse out the necessary precursor data for the radial_distance_assessment class above. Ideally, should operate remotely on an IPython engine to allow for an asychronous parallel workflow with each replicate (universe object) analyzed simultaneously on a different core.'''
    PPCH_PO4_selection = universe_object.select_atoms('resname PPCH and name PO4')
    FORS_AM2_selection = universe_object.select_atoms('resname FORS and name AM2')
    CHOL_ROH_selection = universe_object.select_atoms('resname CHOL and name ROH')
    remaining_headgroup_selection = universe_object.select_atoms('(resname DOPE or resname DOPX or resname POPS) and name PO4')
    total_PPCH_PO4_particles = PPCH_PO4_selection.n_atoms
    total_FORS_AM2_particles = FORS_AM2_selection.n_atoms
    total_CHOL_ROH_particles = CHOL_ROH_selection.n_atoms
    total_remaining_particles = remaining_headgroup_selection.n_atoms
    if FORS_present:
        threshold = 275 #smaller threshold in presence of FORS virions (which are smaller)
    elif control_condition:
        threshold = 600 #60 nm threshold for control conditions
    else:
        threshold = 285 #28.5 nm threshold for outer leaflet PPCH PO4 (currently testing)
    all_lipid_selection = universe_object.select_atoms('resname PPCH or resname CHOL or resname POPS or resname DOPX or resname DOPE or resname FORS')

    list_min_PPCH_PO4_distances = []
    list_max_PPCH_PO4_distances = []
    list_average_PPCH_PO4_distances = []
    list_std_dev_PPCH_PO4_distances = []
    list_frame_numbers = []
    list_PPCH_percent_above_treshold = []

    list_min_FORS_AM2_distances = []
    list_max_FORS_AM2_distances = []
    list_average_FORS_AM2_distances = []
    list_std_dev_FORS_AM2_distances = []
    list_frame_numbers = []
    list_FORS_percent_above_treshold = []

    list_min_CHOL_ROH_distances = []
    list_max_CHOL_ROH_distances = []
    list_average_CHOL_ROH_distances = []
    list_std_dev_CHOL_ROH_distances = []
    list_CHOL_ROH_midpoint_distances = []
    list_CHOL_ROH_percent_above_threshold = []
    list_CHOL_ROH_percent_below_threshold = []

    list_min_remaining_headgroup_distances = []
    list_max_remaining_headgroup_distances = []
    list_average_remaining_headgroup_distances = []
    list_std_dev_remaining_headgroup_distances = []
    list_remaining_headgroup_midpoint_distances = []
    list_remaining_headgroup_percent_above_threshold = []
    list_remaining_headgroup_percent_below_threshold = []

    for ts in universe_object.trajectory[::100]: #every 100th frame (roughly every 10 ns)
        PPCH_PO4_coordinates = PPCH_PO4_selection.coordinates()
        if FORS_present:
            FORS_AM2_coordinates = FORS_AM2_selection.coordinates()
        CHOL_ROH_coordinates = CHOL_ROH_selection.coordinates()
        remaining_headgroup_coordinates = remaining_headgroup_selection.coordinates()
        all_lipid_centroid = all_lipid_selection.centroid()
        PPCH_PO4_coordinates -= all_lipid_centroid #place the centroid of the vesicle at the origin
        if FORS_present:
            FORS_AM2_coordinates -= all_lipid_centroid
        CHOL_ROH_coordinates -= all_lipid_centroid #place the centroid of the vesicle at the origin
        remaining_headgroup_coordinates -= all_lipid_centroid
        spherical_polar_PPCH_PO4_coordinates = convert_cartesian_array_to_spherical_array(PPCH_PO4_coordinates)
        if FORS_present:
            spherical_polar_FORS_AM2_coordinates = convert_cartesian_array_to_spherical_array(FORS_AM2_coordinates)
        spherical_polar_CHOL_ROH_coordinates = convert_cartesian_array_to_spherical_array(CHOL_ROH_coordinates)
        spherical_polar_remaining_coordinates = convert_cartesian_array_to_spherical_array(remaining_headgroup_coordinates)
        #determine the minimum and maximum and average radial distance so I get an idea of leaflet distribution for PPCH PO4
        minimum_PPCH_PO4_radial_distance = spherical_polar_PPCH_PO4_coordinates[...,0].min()
        maximum_PPCH_PO4_radial_distance = spherical_polar_PPCH_PO4_coordinates[...,0].max()
        average_PPCH_PO4_radial_distance = numpy.average(spherical_polar_PPCH_PO4_coordinates[...,0])
        std_dev_PPCH_PO4_radial_distance = numpy.std(spherical_polar_PPCH_PO4_coordinates[...,0])
        number_of_PPCH_PO4_radial_distances_above_threshold = (spherical_polar_PPCH_PO4_coordinates[...,0] > threshold).sum()
        percent_PPCH_PO4_above_treshold = float(number_of_PPCH_PO4_radial_distances_above_threshold) / float(total_PPCH_PO4_particles) * 100.
        list_PPCH_percent_above_treshold.append(percent_PPCH_PO4_above_treshold)
        list_min_PPCH_PO4_distances.append(minimum_PPCH_PO4_radial_distance)
        list_max_PPCH_PO4_distances.append(maximum_PPCH_PO4_radial_distance)
        list_average_PPCH_PO4_distances.append(average_PPCH_PO4_radial_distance)
        list_std_dev_PPCH_PO4_distances.append(std_dev_PPCH_PO4_radial_distance)
        #since I'm currently applying the same cutoff threshold for FORS in the analysis proper, treat FORS AM2 the same way as PPCH PO4 (I think this may actually be problematic, but check the results first)
        #determine the minimum and maximum and average radial distance so I get an idea of leaflet distribution for PPCH PO4
        if FORS_present:
            minimum_FORS_AM2_radial_distance = spherical_polar_FORS_AM2_coordinates[...,0].min()
            maximum_FORS_AM2_radial_distance = spherical_polar_FORS_AM2_coordinates[...,0].max()
            average_FORS_AM2_radial_distance = numpy.average(spherical_polar_FORS_AM2_coordinates[...,0])
            std_dev_FORS_AM2_radial_distance = numpy.std(spherical_polar_FORS_AM2_coordinates[...,0])
            number_of_FORS_AM2_radial_distances_above_threshold = (spherical_polar_FORS_AM2_coordinates[...,0] > threshold).sum()
            percent_FORS_AM2_above_treshold = float(number_of_FORS_AM2_radial_distances_above_threshold) / float(total_FORS_AM2_particles) * 100.
            list_FORS_percent_above_treshold.append(percent_FORS_AM2_above_treshold)
            list_min_FORS_AM2_distances.append(minimum_FORS_AM2_radial_distance)
            list_max_FORS_AM2_distances.append(maximum_FORS_AM2_radial_distance)
            list_average_FORS_AM2_distances.append(average_FORS_AM2_radial_distance)
            list_std_dev_FORS_AM2_distances.append(std_dev_FORS_AM2_radial_distance)
        #similarly for CHOL ROH
        minimum_CHOL_ROH_radial_distance = spherical_polar_CHOL_ROH_coordinates[...,0].min()
        #maximum_CHOL_ROH_radial_distance = spherical_polar_CHOL_ROH_coordinates[...,0].max()
        maximum_CHOL_ROH_radial_distance = numpy.sort(spherical_polar_CHOL_ROH_coordinates[...,0])[-2] #it seems we have a CHOL floater for part of the simulation so avoid that one
        midpoint_CHOL_ROH_radial_distance = numpy.average(numpy.array([minimum_CHOL_ROH_radial_distance,maximum_CHOL_ROH_radial_distance])) #unbiased midpoint for CHOL for upper / lower leaflet cutoff
        average_CHOL_ROH_radial_distance = numpy.average(spherical_polar_CHOL_ROH_coordinates[...,0])
        std_dev_CHOL_ROH_radial_distance = numpy.std(spherical_polar_CHOL_ROH_coordinates[...,0])
        list_min_CHOL_ROH_distances.append(minimum_CHOL_ROH_radial_distance)
        list_max_CHOL_ROH_distances.append(maximum_CHOL_ROH_radial_distance)
        list_average_CHOL_ROH_distances.append(average_CHOL_ROH_radial_distance)
        list_std_dev_CHOL_ROH_distances.append(std_dev_CHOL_ROH_radial_distance)
        list_CHOL_ROH_midpoint_distances.append(midpoint_CHOL_ROH_radial_distance)
        #for CHOL the threshold will be the unbiased midpoint (except for control conditions, which have absolute lipid leaflet assignments)
        if not control_condition:
            number_of_CHOL_ROH_radial_distances_above_threshold = (spherical_polar_CHOL_ROH_coordinates[...,0] > midpoint_CHOL_ROH_radial_distance).sum()
            number_of_CHOL_ROH_radial_distances_below_threshold = (spherical_polar_CHOL_ROH_coordinates[...,0] < midpoint_CHOL_ROH_radial_distance).sum()
        else:
            number_of_CHOL_ROH_radial_distances_above_threshold = (spherical_polar_CHOL_ROH_coordinates[...,0] > threshold).sum()
            number_of_CHOL_ROH_radial_distances_below_threshold = (spherical_polar_CHOL_ROH_coordinates[...,0] < threshold).sum()
        percent_CHOL_ROH_above_treshold = float(number_of_CHOL_ROH_radial_distances_above_threshold) / float(total_CHOL_ROH_particles) * 100.
        percent_CHOL_ROH_below_treshold = float(number_of_CHOL_ROH_radial_distances_below_threshold) / float(total_CHOL_ROH_particles) * 100.
        list_CHOL_ROH_percent_above_threshold.append(percent_CHOL_ROH_above_treshold)
        list_CHOL_ROH_percent_below_threshold.append(percent_CHOL_ROH_below_treshold)
        #similarly for the remaining headgroup particles (for DOPE/X and POPS)
        minimum_remaining_headgroup_radial_distance = spherical_polar_remaining_coordinates[...,0].min()
        maximum_remaining_headgroup_radial_distance = numpy.sort(spherical_polar_remaining_coordinates[...,0])[-2] #dealing with floater(s)
        midpoint_remaining_headgroup_radial_distance = numpy.average(numpy.array([minimum_remaining_headgroup_radial_distance,maximum_remaining_headgroup_radial_distance])) #unbiased midpoint for upper / lower leaflet cutoff
        average_remaining_headgroup_radial_distance = numpy.average(spherical_polar_remaining_coordinates[...,0])
        std_dev_remaining_headgroup_radial_distance = numpy.std(spherical_polar_remaining_coordinates[...,0])

        if not control_condition:
            number_of_remaining_headgroup_radial_distances_above_threshold = (spherical_polar_remaining_coordinates[...,0] > midpoint_remaining_headgroup_radial_distance).sum()
            number_of_remaining_headgroup_radial_distances_below_threshold = (spherical_polar_remaining_coordinates[...,0] < midpoint_remaining_headgroup_radial_distance).sum()
        else:
            number_of_remaining_headgroup_radial_distances_above_threshold = (spherical_polar_remaining_coordinates[...,0] > threshold).sum()
            number_of_remaining_headgroup_radial_distances_below_threshold = (spherical_polar_remaining_coordinates[...,0] < threshold).sum()
        percent_remaining_headgroup_above_threshold = float(number_of_remaining_headgroup_radial_distances_above_threshold) / float(total_remaining_particles) * 100.
        percent_remaining_headgroup_below_threshold = float(number_of_remaining_headgroup_radial_distances_below_threshold) / float(total_remaining_particles) * 100.
        list_min_remaining_headgroup_distances.append(minimum_remaining_headgroup_radial_distance)
        list_max_remaining_headgroup_distances.append(maximum_remaining_headgroup_radial_distance)
        list_average_remaining_headgroup_distances.append(average_remaining_headgroup_radial_distance)
        list_std_dev_remaining_headgroup_distances.append(std_dev_remaining_headgroup_radial_distance)
        list_remaining_headgroup_midpoint_distances.append(midpoint_remaining_headgroup_radial_distance)
        list_remaining_headgroup_percent_above_threshold.append(percent_remaining_headgroup_above_threshold)
        list_remaining_headgroup_percent_below_threshold.append(percent_remaining_headgroup_below_threshold)
        
        
        frame_number = ts.frame
        list_frame_numbers.append(frame_number)

    if not FORS_present:
        return (list_min_PPCH_PO4_distances,list_max_PPCH_PO4_distances,list_average_PPCH_PO4_distances,list_std_dev_PPCH_PO4_distances,list_frame_numbers,list_PPCH_percent_above_treshold,list_min_CHOL_ROH_distances,list_max_CHOL_ROH_distances,list_average_CHOL_ROH_distances,list_std_dev_CHOL_ROH_distances,list_CHOL_ROH_midpoint_distances,list_CHOL_ROH_percent_above_threshold,list_CHOL_ROH_percent_below_threshold,list_min_remaining_headgroup_distances,list_max_remaining_headgroup_distances,list_average_remaining_headgroup_distances,list_std_dev_remaining_headgroup_distances,list_remaining_headgroup_midpoint_distances,list_remaining_headgroup_percent_above_threshold,list_remaining_headgroup_percent_below_threshold,threshold)
    else:
        return (list_min_FORS_AM2_distances,list_max_FORS_AM2_distances,list_average_FORS_AM2_distances,list_std_dev_FORS_AM2_distances,list_frame_numbers,list_FORS_percent_above_treshold,list_min_PPCH_PO4_distances,list_max_PPCH_PO4_distances,list_average_PPCH_PO4_distances,list_std_dev_PPCH_PO4_distances,list_PPCH_percent_above_treshold,list_min_CHOL_ROH_distances,list_max_CHOL_ROH_distances,list_average_CHOL_ROH_distances,list_std_dev_CHOL_ROH_distances,list_CHOL_ROH_midpoint_distances,list_CHOL_ROH_percent_above_threshold,list_CHOL_ROH_percent_below_threshold,list_min_remaining_headgroup_distances,list_max_remaining_headgroup_distances,list_average_remaining_headgroup_distances,list_std_dev_remaining_headgroup_distances,list_remaining_headgroup_midpoint_distances,list_remaining_headgroup_percent_above_threshold,list_remaining_headgroup_percent_below_threshold,threshold)
        
def create_control_universe_data(flu_coordinate_file_path, output_path='/sansom/n22/bioc1009/spherical_Voronoi_virus_work/', num_frames=5000):
    '''Take a flu simulation coordinate file as input and output two different control xtc files for my ipynb Voronoi analysis (area per lipid) workflow. Will probably set the two control xtc files to have lipid densities (and therefore average area per lipid values) that differ by a factor of two. Will assume the input flu simulation data does NOT include FORS.'''
    import MDAnalysis.coordinates.XTC
    import scipy
    import scipy.spatial.distance
    input_flu_coordinate_file_universe_object = MDAnalysis.Universe(flu_coordinate_file_path)
    dict_lipid_residue_data = {'DOPX':{'sel_string':'resname DOPX and name PO4'},'DOPE':{'sel_string':'resname DOPE and name PO4'},'POPS':{'sel_string':'resname POPS and name PO4'},'CHOL':{'sel_string':'resname CHOL and name ROH'},'PPCH':{'sel_string':'resname PPCH and name PO4'}}
    total_residue_headgroup_coordinates_inner_leaflet = 0
    total_residue_headgroup_coordinates_outer_leaflet = 0
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        sel_string = residue_subdictionary['sel_string']
        selection = input_flu_coordinate_file_universe_object.select_atoms(sel_string)
        residue_subdictionary['selection'] = selection
        if residue_name in ['DOPX','DOPE','POPS']:
            total_residue_headgroup_coordinates_inner_leaflet += residue_subdictionary['selection'].n_atoms
        else: 
            total_residue_headgroup_coordinates_outer_leaflet += residue_subdictionary['selection'].n_atoms
    #now, I want to generate a pseudo random distribution of points on two spheres (leaflets) of different radii -- the number of points should match up nicely with the number of residue headgroup coords above [though I will make a second data set that removes about half the points I think]
    prng = numpy.random.RandomState(117) 
    inner_radius = 500
    outer_radius = 700
    inner_leaflet_coord_array = generate_random_array_spherical_generators(total_residue_headgroup_coordinates_inner_leaflet,inner_radius,prng)
    outer_leaflet_coord_array = generate_random_array_spherical_generators(total_residue_headgroup_coordinates_outer_leaflet,outer_radius,prng)
    #ensure that none of the points are pathologically close
    inner_dist_array = scipy.spatial.distance.pdist(inner_leaflet_coord_array)
    outer_dist_array = scipy.spatial.distance.pdist(outer_leaflet_coord_array)
    assert inner_dist_array.min() > 0.01, "Random sphere generators are pathologically close."
    assert outer_dist_array.min() > 0.01, "Random sphere generators are pathologically close."
    #split up the inner and outer random coordinate arrays and set the corresponding AtomGroup coords to these random positions
    inner_leaflet_particle_counter = 0
    outer_leaflet_particle_counter = 0
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        if residue_name in ['DOPX','DOPE','POPS']:
            num_atoms = residue_subdictionary['selection'].n_atoms
            residue_subdictionary['selection'].set_positions(inner_leaflet_coord_array[inner_leaflet_particle_counter:inner_leaflet_particle_counter + num_atoms,...])
            inner_leaflet_particle_counter += num_atoms
        else: #outer leaflet
            num_atoms = residue_subdictionary['selection'].n_atoms
            residue_subdictionary['selection'].set_positions(outer_leaflet_coord_array[outer_leaflet_particle_counter:outer_leaflet_particle_counter + num_atoms,...])
            outer_leaflet_particle_counter += num_atoms
    #now write the first control xtc file with the above random positions on sphere surface
    xtc_writer_instace_1 = MDAnalysis.coordinates.XTC.XTCWriter(output_path + 'control_traj_1.xtc',total_residue_headgroup_coordinates_outer_leaflet + total_residue_headgroup_coordinates_inner_leaflet)
    frames_to_write = num_frames
    while frames_to_write > 0:
        xtc_writer_instace_1.write(input_flu_coordinate_file_universe_object.select_atoms('(resname DOPX and name PO4) or (resname DOPE and name PO4) or (resname POPS and name PO4) or (resname CHOL and name ROH) or (resname PPCH and name PO4)')) #5000 frames with the same custom random coordinates
        frames_to_write -= 1
    #now, set up for writing a second control xtc file, with about half as many total coordinates in each leaflet
    merged_halved_atom_groups = None
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        full_selection_current_residue_type = residue_subdictionary['selection']
        num_atoms_current_residue_type = full_selection_current_residue_type.n_atoms
        approx_half_num_atoms_current_residue_type = int(float(num_atoms_current_residue_type)/2.)
        halved_atomgroup_current_residue_type = full_selection_current_residue_type[0:approx_half_num_atoms_current_residue_type]
        if not merged_halved_atom_groups:
            merged_halved_atom_groups = halved_atomgroup_current_residue_type
        else: #start concatenating once initialized
            merged_halved_atom_groups += halved_atomgroup_current_residue_type
    #now write the second control xtc file with approx. half as many coordinates in each leaflet
    xtc_writer_instace_2 = MDAnalysis.coordinates.XTC.XTCWriter(output_path + 'control_traj_2.xtc',merged_halved_atom_groups.n_atoms)
    frames_to_write = num_frames
    while frames_to_write > 0:
        xtc_writer_instace_2.write(merged_halved_atom_groups) #5000 frames with the same custom random coordinates
        frames_to_write -= 1



def create_control_universe_coord_data(flu_coordinate_file_path,output_path='/sansom/n22/bioc1009/spherical_Voronoi_virus_work/'):
    '''Produce the .gro coordinate files that correspond to the .xtc files produced by the similarly-named function. This should eventually be merged into the matching xtc producer function, but writing this so I don't have to rerun that slow code for now'''
    import MDAnalysis.coordinates.GRO
    import scipy
    import scipy.spatial.distance
    input_flu_coordinate_file_universe_object = MDAnalysis.Universe(flu_coordinate_file_path)
    dict_lipid_residue_data = {'DOPX':{'sel_string':'resname DOPX and name PO4'},'DOPE':{'sel_string':'resname DOPE and name PO4'},'POPS':{'sel_string':'resname POPS and name PO4'},'CHOL':{'sel_string':'resname CHOL and name ROH'},'PPCH':{'sel_string':'resname PPCH and name PO4'}}
    total_residue_headgroup_coordinates_inner_leaflet = 0
    total_residue_headgroup_coordinates_outer_leaflet = 0
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        sel_string = residue_subdictionary['sel_string']
        selection = input_flu_coordinate_file_universe_object.select_atoms(sel_string)
        residue_subdictionary['selection'] = selection
        if residue_name in ['DOPX','DOPE','POPS']:
            total_residue_headgroup_coordinates_inner_leaflet += residue_subdictionary['selection'].n_atoms
        else: 
            total_residue_headgroup_coordinates_outer_leaflet += residue_subdictionary['selection'].n_atoms
    #now, I want to generate a pseudo random distribution of points on two spheres (leaflets) of different radii -- the number of points should match up nicely with the number of residue headgroup coords above [though I will make a second data set that removes about half the points I think]
    prng = numpy.random.RandomState(117) 
    inner_radius = 500
    outer_radius = 700
    inner_leaflet_coord_array = generate_random_array_spherical_generators(total_residue_headgroup_coordinates_inner_leaflet,inner_radius,prng)
    outer_leaflet_coord_array = generate_random_array_spherical_generators(total_residue_headgroup_coordinates_outer_leaflet,outer_radius,prng)
    #ensure that none of the points are pathologically close
    inner_dist_array = scipy.spatial.distance.pdist(inner_leaflet_coord_array)
    outer_dist_array = scipy.spatial.distance.pdist(outer_leaflet_coord_array)
    assert inner_dist_array.min() > 0.01, "Random sphere generators are pathologically close."
    assert outer_dist_array.min() > 0.01, "Random sphere generators are pathologically close."
    #split up the inner and outer random coordinate arrays and set the corresponding AtomGroup coords to these random positions
    inner_leaflet_particle_counter = 0
    outer_leaflet_particle_counter = 0
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        if residue_name in ['DOPX','DOPE','POPS']:
            num_atoms = residue_subdictionary['selection'].n_atoms
            residue_subdictionary['selection'].set_positions(inner_leaflet_coord_array[inner_leaflet_particle_counter:inner_leaflet_particle_counter + num_atoms,...])
            inner_leaflet_particle_counter += num_atoms
        else: #outer leaflet
            num_atoms = residue_subdictionary['selection'].n_atoms
            residue_subdictionary['selection'].set_positions(outer_leaflet_coord_array[outer_leaflet_particle_counter:outer_leaflet_particle_counter + num_atoms,...])
            outer_leaflet_particle_counter += num_atoms
    #now write the first control gro file with the above random positions on sphere surface
    gro_writer_instace_1 = MDAnalysis.coordinates.GRO.GROWriter(output_path + 'control_traj_1.gro')
    gro_writer_instace_1.write(input_flu_coordinate_file_universe_object.select_atoms('(resname DOPX and name PO4) or (resname DOPE and name PO4) or (resname POPS and name PO4) or (resname CHOL and name ROH) or (resname PPCH and name PO4)')) 
    #now, set up for writing a second control gro file, with about half as many total coordinates in each leaflet
    merged_halved_atom_groups = None
    for residue_name, residue_subdictionary in dict_lipid_residue_data.iteritems():
        full_selection_current_residue_type = residue_subdictionary['selection']
        num_atoms_current_residue_type = full_selection_current_residue_type.n_atoms
        approx_half_num_atoms_current_residue_type = int(float(num_atoms_current_residue_type)/2.)
        halved_atomgroup_current_residue_type = full_selection_current_residue_type[0:approx_half_num_atoms_current_residue_type]
        if not merged_halved_atom_groups:
            merged_halved_atom_groups = halved_atomgroup_current_residue_type
        else: #start concatenating once initialized
            merged_halved_atom_groups += halved_atomgroup_current_residue_type
    #now write the second control gro file with approx. half as many coordinates in each leaflet
    gro_writer_instace_2 = MDAnalysis.coordinates.GRO.GROWriter(output_path + 'control_traj_2.gro')
    gro_writer_instace_2.write(merged_halved_atom_groups) 

def TMD_particle_selector_dengue(asymmetric_unit_input_array):
    '''The asymmetric_unit_input_array should be a single numpy array of coordinates for a given asymmetric unit of the dengue viron. This function assumes an asymmetric unit topology of 3 x E followed by 3 x M. Returns the centroids of the individual TMDs: [list_E_TMD_centroids, list_M_TMD_centroids].'''
    E_TMD_1_coordinate_array = asymmetric_unit_input_array[977:1013,...]
    E_TMD_2_coordinate_array = asymmetric_unit_input_array[2047:2083,...]
    E_TMD_3_coordinate_array = asymmetric_unit_input_array[3117:3153,...]
    E_TMD_4_coordinate_array = asymmetric_unit_input_array[1027:1061,...]
    E_TMD_5_coordinate_array = asymmetric_unit_input_array[2097:2131,...]
    E_TMD_6_coordinate_array = asymmetric_unit_input_array[3167:3201,...]
    M_TMD_1_coordinate_array = asymmetric_unit_input_array[3299:3327,...]
    M_TMD_2_coordinate_array = asymmetric_unit_input_array[3457:3485,...]
    M_TMD_3_coordinate_array = asymmetric_unit_input_array[3615:3643,...]
    M_TMD_4_coordinate_array = asymmetric_unit_input_array[3336:3363,...]
    M_TMD_5_coordinate_array = asymmetric_unit_input_array[3494:3521,...]
    M_TMD_6_coordinate_array = asymmetric_unit_input_array[3652:3679,...]
    list_E_TMD_centroids = [numpy.average(E_TMD_array,axis=0) for E_TMD_array in [E_TMD_1_coordinate_array,E_TMD_2_coordinate_array,E_TMD_3_coordinate_array,E_TMD_4_coordinate_array,E_TMD_5_coordinate_array,E_TMD_6_coordinate_array]]
    list_M_TMD_centroids = [numpy.average(M_TMD_array,axis=0) for M_TMD_array in [M_TMD_1_coordinate_array,M_TMD_2_coordinate_array,M_TMD_3_coordinate_array,M_TMD_4_coordinate_array,M_TMD_5_coordinate_array,M_TMD_6_coordinate_array]]
    return [list_E_TMD_centroids,list_M_TMD_centroids]

