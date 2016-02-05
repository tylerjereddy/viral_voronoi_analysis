import unittest
import voronoi_utility
import numpy as np
import scipy
from scipy.spatial import SphericalVoronoi
import numpy.testing
import voronoi_analysis_library
import MDAnalysis
import MDAnalysis.coordinates
from MDAnalysis.coordinates.xdrfile.XTC import XTCWriter
from testfixtures import TempDirectory

class TestRawNeighbourAnalysis(unittest.TestCase):

    def setUp(self):
        self.small_data_dict = {'POPC':{},'PPCE':{}} #build a dict with simple coords to test neighbour properties
        self.template_array = np.zeros((3,3))
        self.a = np.copy(self.template_array)
        self.b = np.copy(self.template_array)
        self.c = np.copy(self.template_array)
        self.a[0,...] = np.ones((1,3))
        self.b[1,...] = np.ones((1,3))
        self.c[2,...] = np.ones((1,3))
        

    def tearDown(self):
        del self.small_data_dict
        del self.template_array
        del self.a
        del self.b
        del self.c

    def test_neighbour_counting_identical_shuffled_coords(self):
        '''a,b,c contain the same vertices in different order and so all Voronoi cells 
        in the system should be counted as having zero neighbours (i.e., the code should
        only count vertex matches as neighbours when < all vertices match).'''
        for resname, subdict in self.small_data_dict.iteritems():
            subdict['voronoi_cell_list_vertex_arrays'] = [[self.a,self.b,self.c]]
            subdict['voronoi_cell_list_vertex_arrays_inner_leaflet'] = [[self.a,self.b,self.c]]
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays')
        inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict = voronoi_neighbour_instance.identify_voronoi_neighbours(frame_index = 0)
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'].keys(), [0], "Inner leaflet POPC should have zero neighbours.") 
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE'].keys(), [0], "Inner leaflet PPCE should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['POPC'].keys(), [0], "Outer leaflet POPC should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'].keys(), [0], "Outer leaflet PPCE should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'][0]['frequency'], 3, "Outer leaflet PPCE should have 3 cells (a,b,c) with 0 neighbours.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'][0]['frequency'], 3, "Inner leaflet POPC should have 3 cells (a,b,c) with 0 neighbours.")

    def test_neighbour_counting_single_unique_cell(self):
        '''Using artificial polygon (a) as a single unique input, with all other artificial polygons containing only 
        one matching vertex to a, a neighbour counting result should indicate that (a) has 4 neighbours (2 from each residue type) and that each of those 4 neighbours has only
        2 neighbours each (1 copy of (a) from each residue), given two residues with three Voronoi cells each.'''
        a = np.array([[3,3,3],[5,5,5],[7,7,7]])
        b = np.array([[3,3,3],[0,0,0],[9,9,9]])
        for resname, subdict in self.small_data_dict.iteritems():
            subdict['voronoi_cell_list_vertex_arrays'] = [[a,b,b]]
            subdict['voronoi_cell_list_vertex_arrays_inner_leaflet'] = [[a,b,b]]
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays')
        inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict = voronoi_neighbour_instance.identify_voronoi_neighbours(frame_index = 0)
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'].keys(), [2,4], "Inner leaflet POPC should have either 2 or 4 neighbours, but got {result}.".format(result=inner_leaflet_neighbour_dict['POPC'].keys())) 
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE'].keys(), [2,4], "Inner leaflet PPCE should have either 2 or 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['POPC'].keys(), [2,4], "Outer leaflet POPC should have either 2 or 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'].keys(), [2,4], "Outer leaflet PPCE should have either 2 or 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'][2]['frequency'], 2, "Outer leaflet PPCE should have 2 cells with 2 neighbours.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'][2]['frequency'], 2, "Inner leaflet POPC should have 2 cells with 2 neighbours.")
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'][4]['frequency'], 1, "Outer leaflet PPCE should have 1 cell with 4 neighbours.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'][4]['frequency'], 1, "Inner leaflet POPC should have 1 cell with 4 neighbours.")

    def test_neighbour_counting_all_unique_cells(self):
        '''Using 3 unique artifical polygons that all share one vertex with each other as input, 
        two residue types with three molecules each should have a neighbour result with
        4 neighbours for all polygons, and frequency of 3 per residue type.'''
        a = np.array([[3,3,3],[5,5,5],[7,7,7]])
        b = np.array([[3,3,3],[0,0,0],[9,9,9]])
        c = np.array([[9,9,9],[1,1,1],[5,5,5]])
        for resname, subdict in self.small_data_dict.iteritems():
            subdict['voronoi_cell_list_vertex_arrays'] = [[a,b,c]]
            subdict['voronoi_cell_list_vertex_arrays_inner_leaflet'] = [[a,b,c]]
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays')
        inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict = voronoi_neighbour_instance.identify_voronoi_neighbours(frame_index = 0)
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'].keys(), [4], "Inner leaflet POPC should have 4 neighbours")
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE'].keys(), [4], "Inner leaflet PPCE should have 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['POPC'].keys(), [4], "Outer leaflet POPC should have 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'].keys(), [4], "Outer leaflet PPCE should have 4 neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'][4]['frequency'], 3, "Outer leaflet PPCE should have 3 cells with 4 neighbours.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'][4]['frequency'], 3, "Inner leaflet POPC should have 3 cells with 4 neighbours.")


class TestSpeciesSpecificNeighbourAnalysis(unittest.TestCase):

    def setUp(self):
        self.small_data_dict = {'POPC':{},'PPCE':{}} #build a dict with simple coords to test neighbour properties
        self.a = np.array([[3,3,3],[5,5,5],[7,7,7]])
        self.b = np.array([[3,3,3],[0,0,0],[9,9,9]])
        self.c = np.array([[9,9,9],[1,1,1],[5,5,5]])

    def tearDown(self):
        del self.small_data_dict
        del self.a
        del self.b
        del self.c

    def test_neighbours_3_unique_polygons(self):
        '''Given 3 unique polygons a,b,c specifying coordinates for Voronoi cells for two residue types,
        and given that each polygon only shares one vertex with any other, the small POPC / POPE system
        should have each polygon in a given residue type paired with two neighbours of the same type and 
        two neighbours of the other type (self should be excluded regardless of residue type).'''
        for resname, subdict in self.small_data_dict.iteritems():
            subdict['voronoi_cell_list_vertex_arrays'] = [[self.a,self.b,self.c]]
            subdict['voronoi_cell_list_vertex_arrays_inner_leaflet'] = [[self.a,self.b,self.c]]
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis_by_type(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays')
        inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict = voronoi_neighbour_instance.identify_voronoi_neighbours(frame_index = 0)
        self.assertEqual(inner_leaflet_neighbour_dict['POPC']['POPC'].keys(), [2], "Inner leaflet POPC should have 2 neighbours of type POPC.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC']['PPCE'].keys(), [2], "Inner leaflet POPC should have 2 neighbours of type PPCE.")
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE']['POPC'].keys(), [2], "Inner leaflet PPCE should have 2 neighbours of type POPC.")
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE']['PPCE'].keys(), [2], "Inner leaflet PPCE should have 2 neighbours of type PPCE.")
        self.assertEqual(outer_leaflet_neighbour_dict['POPC']['POPC'].keys(), [2], "Outer leaflet POPC should have 2 neighbours of type POPC.")
        self.assertEqual(outer_leaflet_neighbour_dict['POPC']['PPCE'].keys(), [2], "Outer leaflet POPC should have 2 neighbours of type PPCE.")
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE']['POPC'].keys(), [2], "Outer leaflet PPCE should have 2 neighbours of type POPC.")
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE']['PPCE'].keys(), [2], "Outer leaflet PPCE should have 2 neighbours of type PPCE.")
        #there should be 3 cells (a,b,c) with two neighbours in all combinations:
        self.assertEqual(len(outer_leaflet_neighbour_dict['PPCE']['PPCE'][2]), 3)
        self.assertEqual(len(inner_leaflet_neighbour_dict['POPC']['PPCE'][2]), 3)


class TestVoronoiAreaSummation(unittest.TestCase):

    def setUp(self):
        self.test_dict = {0: 1.55, 1: 1.32, 2: 1.56, 3: 0.97, 4: 0.55, 5: 2.9}
        #arbitrarily set the first half to protein cells
        self.total_protein_area = 1.55 + 1.32 + 1.56
        self.total_lipid_area = 0.97 + 0.55 + 2.9

    def tearDown(self):
        del self.test_dict
        del self.total_protein_area
        del self.total_lipid_area

    def test_summation(self):
        '''Confirm proper summation of protein and lipid surface areas in sum_Voronoi_cell_surface_areas() function.'''
        protein_SA_sum, lipid_SA_sum = voronoi_analysis_library.sum_Voronoi_cell_surface_areas(0,3,self.test_dict)
        self.assertEqual(protein_SA_sum, self.total_protein_area, "Sum of protein surface areas incorrect.")
        self.assertEqual(lipid_SA_sum, self.total_lipid_area, "Sum of lipid surface areas incorrect.")


class TestDengueTMDSelector(unittest.TestCase):

    def setUp(self):
        self.test_array = np.random.random_sample((4000,3))
        self.third_E_TMD_centroid = np.average(self.test_array[3117:3153,...], axis = 0)
        self.sixth_M_TMD_centroid = np.average(self.test_array[3652:3679,...], axis = 0)
        
    def tearDown(self):
        del self.test_array
        del self.third_E_TMD_centroid
        del self.sixth_M_TMD_centroid

    def test_selector_return_lists(self):
        list_E_TMD_centroids, list_M_TMD_centroids = voronoi_analysis_library.TMD_particle_selector_dengue(self.test_array)
        np.testing.assert_array_almost_equal(list_E_TMD_centroids[2], self.third_E_TMD_centroid, err_msg = "Did not retrieve correct centroid for E TMD #3.")
        np.testing.assert_array_almost_equal(list_M_TMD_centroids[-1], self.sixth_M_TMD_centroid, err_msg = "Did not retrieve correct centroid for M TMD #6.")


class TestFluTMDSelector(unittest.TestCase):

    def setUp(self):
        self.test_array = np.random.random_sample((4000,3))
        
    def tearDown(self):
        del self.test_array

    def test_HA(self):
        output_array = voronoi_analysis_library.TMD_particle_selector(self.test_array, 'HA')
        self.assertEqual(output_array.shape, (3,3), "The array of HA TMD centroid coordinates should have shape (3,3)")
        HA_TMD_2_centroid = np.average(self.test_array[2305:2369],axis = 0)
        np.testing.assert_array_almost_equal(output_array[1,...], HA_TMD_2_centroid, decimal=6, err_msg="The second HA TMD centroid is incorrect.")

    def test_NA(self):
        output_array = voronoi_analysis_library.TMD_particle_selector(self.test_array, 'NA')
        self.assertEqual(output_array.shape, (4,3), "The array of NA TMD centroid coordinates should have shape (4,3)")
        NA_TMD_4_centroid = np.average(self.test_array[3083:3140],axis = 0)
        np.testing.assert_array_almost_equal(output_array[3,...], NA_TMD_4_centroid, decimal=6, err_msg="The fourth NA TMD centroid is incorrect.")

    def test_M2(self):
        output_array = voronoi_analysis_library.TMD_particle_selector(self.test_array, 'M2')
        self.assertEqual(output_array.shape, (4,3), "The array of M2 TMD centroid coordinates should have shape (4,3)")
        M2_TMD_1_centroid = np.average(self.test_array[13:64],axis = 0)
        np.testing.assert_array_almost_equal(output_array[0,...], M2_TMD_1_centroid, decimal=6, err_msg="The first M2 TMD centroid is incorrect.")

    def test_lipid(self):
        output_array = voronoi_analysis_library.TMD_particle_selector(self.test_array, 'lipid')
        np.testing.assert_array_almost_equal(output_array, self.test_array, decimal=6, err_msg="Lipid output array should be unchanged from input.")
        
class TestSurfaceAreaSphere(unittest.TestCase):

    def setUp(self):
        self.sphere_radius = 1.997842
        self.google_result = 50.157 #SA from google calculation

    def tearDown(self):
        del self.sphere_radius
        del self.google_result

    def test_SA_calculation(self):
        calculated_SA = voronoi_analysis_library.calculate_surface_area_sphere(self.sphere_radius)
        self.assertAlmostEqual(calculated_SA, self.google_result, places = 3, msg = "Calculated surface area of sphere is not correct to 3 decimal places.")

class TestVoronoiAnalysisLoop(unittest.TestCase):
    '''Functional test(s) to ensure stability / correctness of the huge voronoi_analysis_library.voronoi_analysis_loop() function.'''

    @classmethod
    def setUpClass(cls):
        cls.u = MDAnalysis.Universe('control_traj_2.gro') #mock flu universe with DOPE/X and POPS inner leaflet; PPCH / CHOL outer leaflet; no protein
        cls.d = TempDirectory()
        #create a short trajectory with the same control coords in each frame
        cls.num_frames = 4
        cls.xtc = cls.d.path + '/control_traj_2_dummy.xtc'
        with XTCWriter(cls.xtc, cls.u.trajectory.n_atoms) as W:
            while cls.num_frames > 0:
                W.write(cls.u)
                cls.num_frames -= 1
        cls.u_multiframe = MDAnalysis.Universe('control_traj_2.gro', cls.xtc)
        cls.loop_result = voronoi_analysis_library.voronoi_analysis_loop(cls.u_multiframe,0,'full',1,control_condition=1)

    @classmethod
    def tearDownClass(cls):
        del cls.u
        cls.d.cleanup()
        del cls.num_frames
        del cls.loop_result
        del cls.xtc
        del cls.u_multiframe

    def test_surface_area_reconstitution(self):
        '''For a control system with no proteins we should achive > 99% surface area reconstitution in all frames.'''
        array_outer_leaflet_percent_SA_reconstitution = np.array(self.loop_result[1])
        array_inner_leaflet_percent_SA_reconstitution = np.array(self.loop_result[2])
        outer_min = array_outer_leaflet_percent_SA_reconstitution.min()
        inner_min = array_inner_leaflet_percent_SA_reconstitution.min()
        self.assertGreaterEqual(outer_min, 99.0, "Outer leaflet % surface area reconsitution drops below 99%. Minimum value found was {mini}.".format(mini=outer_min))
        self.assertGreaterEqual(inner_min, 99.0, "Inner leaflet % surface area reconsitution drops below 99%. Minimum value found was {mini}.".format(mini=inner_min))

    def test_frame_count(self):
        '''Ensure that data structures returned by voronoi_analysis_loop match the number of frames in the input data.'''
        list_frame_numbers = self.loop_result[0]
        list_outer_leaflet_percent_SA_reconstitution = self.loop_result[1]
        list_inner_leaflet_percent_SA_reconstitution = self.loop_result[2]
        self.assertEqual(list_frame_numbers, [0,1,2,3], "Incorrect frame numbers returned by voronoi_analysis_loop: {returned_list}.".format(returned_list=list_frame_numbers))
        tuple_list_lengths = len(list_outer_leaflet_percent_SA_reconstitution), len(list_inner_leaflet_percent_SA_reconstitution)
        self.assertEqual(tuple_list_lengths, (4,4), "Surface area reconstitution data structures are inconsistent with the number of frames in the data received by voronoi_analysis_loop. Actual list lengths received: {list_tuple}".format(list_tuple=tuple_list_lengths))

    def test_dictionary_headgroup_data_structure(self):
        '''Check data structures in the dictionary_headgroup_data as returned by voronoi_analysis_loop.'''
        dictionary_headgroup_data = self.loop_result[3]
        keys = dictionary_headgroup_data.keys()
        self.assertEqual(keys, ['POPS', 'DOPE', 'CHOL', 'PPCH', 'DOPX'], "Incorrect set of keys in dictionary_headgroup_data: {keys}".format(keys=keys))
        self.assertEqual(len(dictionary_headgroup_data['CHOL']['voronoi_cell_avg_values_list']),4,"List of outer leaflet Voronoi cell areas is not consistent with the number of frames in the trajectory input to voronoi_analysis_loop.")

    def test_voronoi_analysis_loop_invariance(self):
        '''As the input test xtc has the same 4 frames, the data extracted from each frame should match (loop invariance).'''
        list_outer_leaflet_percent_SA_reconstitution = self.loop_result[1]
        list_inner_leaflet_percent_SA_reconstitution = self.loop_result[2]
        self.assertTrue(list_outer_leaflet_percent_SA_reconstitution.count(list_outer_leaflet_percent_SA_reconstitution[0]) == len(list_outer_leaflet_percent_SA_reconstitution), "Not all outer leaflet surface area reconstitution values are the same for each frame, despite identical coordinates for each input frame in test trajectory.")
        self.assertTrue(list_inner_leaflet_percent_SA_reconstitution.count(list_inner_leaflet_percent_SA_reconstitution[0]) == len(list_inner_leaflet_percent_SA_reconstitution), "Not all inner leaflet surface area reconstitution values are the same for each frame, despite identical coordinates for each input frame in test trajectory.")

class TestVoronoiAreaDict(unittest.TestCase):

    def setUp(self):
        self.prng = numpy.random.RandomState(117) 
        self.sphere_radius = 17.3
        self.num_points = 30
        self.random_coords_sphere = voronoi_utility.generate_random_array_spherical_generators(self.num_points, self.sphere_radius, self.prng)

    def tearDown(self):
        del self.prng
        del self.sphere_radius
        del self.num_points 
        del self.random_coords_sphere 

    def test_voronoi_dict_data_struct(self):
        '''Test the basic data structure returned by produce_Voronoi_area_dict() function.'''
        voronoi_instance = SphericalVoronoi(self.random_coords_sphere,self.sphere_radius)
        voronoi_instance.sort_vertices_of_regions()
        list_voronoi_polygon_vertex_indices = voronoi_instance.regions 
        list_voronoi_polygon_vertices = []
        for region in list_voronoi_polygon_vertex_indices:
            list_voronoi_polygon_vertices.append(voronoi_instance.vertices[region])
        dictionary_voronoi_polygon_surface_areas = voronoi_analysis_library.produce_Voronoi_area_dict(list_voronoi_polygon_vertices,self.sphere_radius)
        self.assertEqual(len(dictionary_voronoi_polygon_surface_areas), self.num_points, "The number of items in dictionary_voronoi_polygon_surface_areas does not match the number of generators in the input test data.")
        final_Voronoi_region_surface_area = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(list_voronoi_polygon_vertices[-1], self.sphere_radius)
        self.assertEqual(dictionary_voronoi_polygon_surface_areas[self.num_points - 1], final_Voronoi_region_surface_area, "The surface area of the final Voronoi region has not been stored correctly in dictionary_voronoi_polygon_surface_areas.")


