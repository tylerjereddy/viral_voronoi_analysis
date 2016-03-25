import unittest
import math
import numpy as np
import scipy
from scipy.spatial import SphericalVoronoi
import numpy.testing
import voronoi_analysis_library
import MDAnalysis
import MDAnalysis.coordinates
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.tests.datafiles import GRO, XTC
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
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays') #using dummy values for the radii as the points aren't actually on a sphere for the test
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
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays') #using dummy values for the radii as the points aren't actually on a sphere for the test
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
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays') #using dummy values for the radii as the points aren't actually on a sphere for the test
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
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis_by_type(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays') #using dummy values for the radii as the points aren't actually on a sphere for the test
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

class CommonTestsVoronoiAnalysisLoop(object):
    '''Tests that are reused for the voronoi_analysis_loop testing.'''

    @classmethod
    def tearDownClass(cls):
        del cls.u
        cls.d.cleanup()
        del cls.num_frames
        del cls.loop_result
        del cls.xtc
        del cls.u_multiframe
        del cls.expected_keys
        del cls.dict_test_key

    def test_frame_count(self):
        '''Ensure that data structures returned by voronoi_analysis_loop match the number of frames in the input data.'''
        list_frame_numbers = self.loop_result[0]
        first_list_percent_SA_reconstitution = self.loop_result[1]
        second_list_percent_SA_reconstitution = self.loop_result[2]
        self.assertEqual(list_frame_numbers, [0,1,2,3], "Incorrect frame numbers returned by voronoi_analysis_loop: {returned_list}.".format(returned_list=list_frame_numbers))
        tuple_list_lengths = len(first_list_percent_SA_reconstitution), len(second_list_percent_SA_reconstitution)
        self.assertEqual(tuple_list_lengths, (4,4), "Surface area reconstitution data structures are inconsistent with the number of frames in the data received by voronoi_analysis_loop. Actual list lengths received: {list_tuple}".format(list_tuple=tuple_list_lengths))

    def test_dictionary_headgroup_data_structure(self):
        '''Check data structures in the dictionary_headgroup_data as returned by voronoi_analysis_loop.'''
        keys = self.dictionary_headgroup_data.keys()
        self.assertEqual(sorted(keys), sorted(self.expected_keys), "Incorrect set of keys in dictionary_headgroup_data: {keys}".format(keys=keys))
        self.assertEqual(len(self.dictionary_headgroup_data[self.dict_test_key]['voronoi_cell_avg_values_list']),4,"List of Voronoi cell areas is not consistent with the number of frames in the trajectory input to voronoi_analysis_loop.")

    def test_voronoi_analysis_loop_invariance(self):
        '''As the input test xtc has the same 4 frames, the data extracted from each frame should match (loop invariance).'''
        list_percent_SA_reconstitution = self.loop_result[1]
        list_percent_SA_reconstitution_2 = self.loop_result[2]
        self.assertTrue(list_percent_SA_reconstitution.count(list_percent_SA_reconstitution[0]) == len(list_percent_SA_reconstitution), "Not all surface area reconstitution values are the same for each frame, despite identical coordinates for each input frame in test trajectory.")
        self.assertTrue(list_percent_SA_reconstitution_2.count(list_percent_SA_reconstitution_2[0]) == len(list_percent_SA_reconstitution_2), "Not all surface area reconstitution values are the same for each frame, despite identical coordinates for each input frame in test trajectory.")

class TestVoronoiAnalysisLoopControl(CommonTestsVoronoiAnalysisLoop,unittest.TestCase):
    '''Functional test(s) to ensure stability / correctness of the huge voronoi_analysis_library.voronoi_analysis_loop() function for CONTROL condition.'''

    @classmethod
    def setUpClass(cls):
        cls.u = MDAnalysis.Universe('control_traj_2_small.gro.bz2') #mock flu universe with DOPE/X and POPS inner leaflet; PPCH / CHOL outer leaflet; no protein
        cls.d = TempDirectory()
        #create short trajectory for testing purposes
        cls.xtc = cls.d.path + '/control_traj_2_dummy.xtc'
        cls.num_frames = 4
        with XTCWriter(cls.xtc, cls.u.trajectory.n_atoms) as W:
            while cls.num_frames > 0:
                W.write(cls.u)
                cls.num_frames -= 1
        cls.u_multiframe = MDAnalysis.Universe('control_traj_2_small.gro.bz2', cls.xtc)
        cls.loop_result = voronoi_analysis_library.voronoi_analysis_loop(cls.u_multiframe,0,'full',1,control_condition=1)
        cls.expected_keys = ['POPS', 'DOPE', 'CHOL', 'PPCH', 'DOPX']
        cls.dictionary_headgroup_data = cls.loop_result[3]
        cls.dict_test_key = 'CHOL'

    def test_surface_area_reconstitution(self):
        '''For a control system with no proteins we should achive > 99% surface area reconstitution in all frames.'''
        array_outer_leaflet_percent_SA_reconstitution = np.array(self.loop_result[1])
        array_inner_leaflet_percent_SA_reconstitution = np.array(self.loop_result[2])
        outer_min = array_outer_leaflet_percent_SA_reconstitution.min()
        inner_min = array_inner_leaflet_percent_SA_reconstitution.min()
        self.assertGreaterEqual(outer_min, 99.0, "Outer leaflet % surface area reconsitution drops below 99%. Minimum value found was {mini}.".format(mini=outer_min))
        self.assertGreaterEqual(inner_min, 99.0, "Inner leaflet % surface area reconsitution drops below 99%. Minimum value found was {mini}.".format(mini=inner_min))

class TestVoronoiAnalysisLoopFlu(CommonTestsVoronoiAnalysisLoop,unittest.TestCase):
    '''Functional test(s) to ensure stability / correctness of the huge voronoi_analysis_library.voronoi_analysis_loop() function for FLU condition.'''

    @classmethod
    def setUpClass(cls):
        cls.u = MDAnalysis.Universe('sim39_final_snapshot_compact_no_solvent_small_proteins_retained.gro.bz2') #actual flu universe
        cls.d = TempDirectory()
        #create short trajectory for testing purposes
        cls.xtc = cls.d.path +  '/flu_dummy.xtc'
        cls.num_frames = 4
        with XTCWriter(cls.xtc, cls.u.trajectory.n_atoms) as W:
            while cls.num_frames > 0:
                W.write(cls.u)
                cls.num_frames -= 1
        cls.u_multiframe = MDAnalysis.Universe('sim39_final_snapshot_compact_no_solvent_small_proteins_retained.gro.bz2', cls.xtc)
        cls.loop_result = voronoi_analysis_library.voronoi_analysis_loop(cls.u_multiframe,0,'full',1,PPCH_PO4_threshold=275,proteins_present='yes',FORS_present='yes')
        cls.expected_keys = ['POPS', 'DOPE', 'CHOL', 'PPCH', 'DOPX', 'FORS', 'protein']
        cls.dictionary_headgroup_data = cls.loop_result[5]
        cls.dict_test_key = 'CHOL'

    def test_surface_area_reconstitution_flu(self):
        '''Various checks for surface area reconstitution properties of the flu virion condition.'''
        sim39_outer_leaflet_lipid_reconstitutions = np.array(self.loop_result[1])
        sim39_outer_leaflet_protein_reconstitutions = np.array(self.loop_result[2])
        sim39_inner_leaflet_lipid_reconstitutions = np.array(self.loop_result[3])
        sim39_inner_leaflet_protein_reconstitutions = np.array(self.loop_result[4])
        
        np.testing.assert_allclose(sim39_outer_leaflet_lipid_reconstitutions + sim39_outer_leaflet_protein_reconstitutions, np.zeros(sim39_outer_leaflet_lipid_reconstitutions.shape) + 100., rtol=1e-07, err_msg="Total % reconstitution of surface area in outer leaflet is not close to 100.")
        np.testing.assert_allclose(sim39_inner_leaflet_lipid_reconstitutions + sim39_inner_leaflet_protein_reconstitutions, np.zeros(sim39_inner_leaflet_lipid_reconstitutions.shape) + 100., rtol=1e-01, err_msg="Total % reconstitution of surface area in inner leaflet is not close to 100.")

        np.testing.assert_array_less(sim39_outer_leaflet_protein_reconstitutions, sim39_outer_leaflet_lipid_reconstitutions, err_msg="For flu simulations, the % SA reconstitution from protein should always be less than the contribution from lipid (outer leaflet).")
        np.testing.assert_array_less(sim39_inner_leaflet_protein_reconstitutions, sim39_inner_leaflet_lipid_reconstitutions, err_msg="For flu simulations, the % SA reconstitution from protein should always be less than the contribution from lipid (inner leaflet).")

        np.testing.assert_array_less(np.zeros(sim39_outer_leaflet_protein_reconstitutions.shape), sim39_outer_leaflet_protein_reconstitutions, err_msg="The % surface area contribution from protein in flu simulations should always exceed 0 (outer leaflet)")
        np.testing.assert_array_less(np.zeros(sim39_inner_leaflet_protein_reconstitutions.shape), sim39_inner_leaflet_protein_reconstitutions, err_msg="The % surface area contribution from protein in flu simulations should always exceed 0 (inner leaflet)")

class TestVoronoiAnalysisLoopDengue(CommonTestsVoronoiAnalysisLoop,unittest.TestCase):
    '''Functional test(s) to ensure stability / correctness of the huge voronoi_analysis_library.voronoi_analysis_loop() function for DENGUE condition.'''

    @classmethod
    def setUpClass(cls):
        cls.u = MDAnalysis.Universe('dengue_final_snapshot_compact_no_solvent.gro.bz2') #dengue universe
        cls.d = TempDirectory()
        #create short trajectory for testing purposes
        cls.xtc = cls.d.path + '/dengue_dummy.xtc'
        cls.num_frames = 4
        with XTCWriter(cls.xtc, cls.u.trajectory.n_atoms) as W:
            while cls.num_frames > 0:
                W.write(cls.u)
                cls.num_frames -= 1
        cls.u_multiframe = MDAnalysis.Universe('dengue_final_snapshot_compact_no_solvent.gro.bz2', cls.xtc)
        cls.loop_result = voronoi_analysis_library.voronoi_analysis_loop(cls.u_multiframe,0,'full',1,proteins_present='yes',dengue_condition=1)
        cls.expected_keys = ['POPC', 'PPCE', 'DPPE', 'CER', 'DUPC', 'DOPS', 'PPCS', 'protein']
        cls.dictionary_headgroup_data = cls.loop_result[5]
        cls.dict_test_key = 'PPCE'

    def test_surface_area_reconstitution_dengue(self):
        '''Various checks for surface area reconstitution properties of the dengue virion condition.'''
        dengue_outer_leaflet_lipid_reconstitutions = np.array(self.loop_result[1])
        dengue_outer_leaflet_protein_reconstitutions = np.array(self.loop_result[2])
        dengue_inner_leaflet_lipid_reconstitutions = np.array(self.loop_result[3])
        dengue_inner_leaflet_protein_reconstitutions = np.array(self.loop_result[4])

        np.testing.assert_allclose(dengue_outer_leaflet_lipid_reconstitutions + dengue_outer_leaflet_protein_reconstitutions, np.zeros(dengue_outer_leaflet_lipid_reconstitutions.shape) + 100., rtol=1e-07, err_msg="Total % reconstitution of surface area in outer leaflet is not close to 100.")
        np.testing.assert_allclose(dengue_inner_leaflet_lipid_reconstitutions + dengue_inner_leaflet_protein_reconstitutions, np.zeros(dengue_inner_leaflet_lipid_reconstitutions.shape) + 100., rtol=1e-01, err_msg="Total % reconstitution of surface area in inner leaflet is not close to 100.")

        np.testing.assert_array_less(dengue_outer_leaflet_protein_reconstitutions, dengue_outer_leaflet_lipid_reconstitutions, err_msg="For dengue simulations, the % SA reconstitution from protein should always be less than the contribution from lipid (outer leaflet).")
        np.testing.assert_array_less(dengue_inner_leaflet_protein_reconstitutions, dengue_inner_leaflet_lipid_reconstitutions, err_msg="For dengue simulations, the % SA reconstitution from protein should always be less than the contribution from lipid (inner leaflet).")

        np.testing.assert_array_less(np.zeros(dengue_outer_leaflet_protein_reconstitutions.shape), dengue_outer_leaflet_protein_reconstitutions, err_msg="The % surface area contribution from protein in dengue simulations should always exceed 0 (outer leaflet)")
        np.testing.assert_array_less(np.zeros(dengue_inner_leaflet_protein_reconstitutions.shape), dengue_inner_leaflet_protein_reconstitutions, err_msg="The % surface area contribution from protein in dengue simulations should always exceed 0 (inner leaflet)")

class TestVoronoiAreaDict(unittest.TestCase):

    def setUp(self):
        self.prng = numpy.random.RandomState(117) 
        self.sphere_radius = 17.3
        self.num_points = 30
        self.random_coords_sphere = voronoi_analysis_library.generate_random_array_spherical_generators(self.num_points, self.sphere_radius, self.prng)

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

class TestSphericalPolygonSA(unittest.TestCase):
    '''Test(s) for calculation of the surface area of a spherical polygon.'''

    @classmethod
    def setUpClass(cls):
        cls.area_unit_sphere = voronoi_analysis_library.calculate_surface_area_sphere(1.0)
        cls.area_larger_sphere = voronoi_analysis_library.calculate_surface_area_sphere(17.0)

    @classmethod
    def tearDownClass(cls):
        del cls.area_unit_sphere
        del cls.area_larger_sphere

    def test_spherical_triangle_SA_eighth_unit_sphere(self):
        '''Test that a spherical triangle covering 1/8 of a unit sphere surface has the appropriate surface area.'''
        # test would fail for the 1/4 case -- perhaps because 3 points are 'co-linear' / equatorial in a sense
        test_polygon = np.array([[0,0,1],[0,1,0],[-1,0,0]]) 
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 1.0)
        theoretical_SA = self.area_unit_sphere / 8.0
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical triangle covering one eighth of unit sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))
        
    def test_spherical_polygon_SA_quarter_unit_sphere(self):
        '''Test that a spherical polygon covering 1/4 of a unit sphere surface has the appropriate surface area.'''
        test_polygon = np.array([[0,0,1],[1,0,0],[0,1,0],[-1,0,0]])
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 1.0)
        theoretical_SA = self.area_unit_sphere / 4.0
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical polygon covering one quarter of unit sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))

    def test_spherical_polygon_SA_3_eighths_unit_sphere(self):
        '''Test that a spherical polygon covering 3/8 of a unit sphere surface has the appropriate surface area.'''
        test_polygon = np.array([[0,0,1],[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]])
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 1.0)
        theoretical_SA = 3./8. * self.area_unit_sphere 
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical polygon covering 3/8  of unit sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))

    def test_spherical_triangle_SA_eighth_large_sphere(self):
        '''Test that a spherical triangle covering 1/8 of a large sphere surface has the appropriate surface area.'''
        # test would fail for the 1/4 case -- perhaps because 3 points are 'co-linear' / equatorial in a sense
        test_polygon = np.array([[0,0,17],[0,17,0],[-17,0,0]]) 
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 17.0)
        theoretical_SA = self.area_larger_sphere / 8.0
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical triangle covering one eighth of larger sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))
        
    def test_spherical_polygon_SA_quarter_larger_sphere(self):
        '''Test that a spherical polygon covering 1/4 of a larger sphere surface has the appropriate surface area.'''
        test_polygon = np.array([[0,0,17],[17,0,0],[0,17,0],[-17,0,0]])
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 17.0)
        theoretical_SA = self.area_larger_sphere / 4.0
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical polygon covering one quarter of larger sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))

    def test_spherical_polygon_SA_3_eighths_larger_sphere(self):
        '''Test that a spherical polygon covering 3/8 of a larger sphere surface has the appropriate surface area.'''
        test_polygon = np.array([[0,0,17],[17,0,0],[0,17,0],[-17,0,0],[0,-17,0]])
        calculated_SA = voronoi_analysis_library.calculate_surface_area_of_a_spherical_Voronoi_polygon(test_polygon, 17.0)
        theoretical_SA = 3./8. * self.area_larger_sphere 
        self.assertAlmostEqual(calculated_SA, theoretical_SA, places = 6, msg="Surface area of spherical polygon covering 3/8  of larger sphere not calculated accurately enough. Calculated: {calc}; Target: {target}".format(calc=calculated_SA, target=theoretical_SA))

class TestHaversineDist(unittest.TestCase):
    '''Unit test(s) for haversine distance calculation on a sphere.'''

    @classmethod
    def setUpClass(cls):
        cls.small_radius = 1.0
        cls.large_radius = 29.892
        cls.small_circumference = 2 * math.pi * cls.small_radius
        cls.large_circumference = 2 * math.pi * cls.large_radius

    @classmethod
    def tearDownClass(cls):
        del cls.small_radius 
        del cls.large_radius
        del cls.small_circumference
        del cls.large_circumference 

    def test_haversine_antipodes_unit_sphere(self):
        '''Test haversine distance calculation for antipodal points on unit sphere.'''
        first_point = np.array([0,0,1])
        second_point = np.array([0,0,-1])
        calculated_distance = voronoi_analysis_library.calculate_haversine_distance_between_spherical_points(first_point, second_point, self.small_radius)
        target_distance = self.small_circumference / 2.
        self.assertAlmostEqual(calculated_distance, target_distance, places=7, msg="The haversine distance between antipodal points on the unit sphere was not calculated correctly. Target: {target}; calculated: {calc}.".format(target=target_distance, calc=calculated_distance))
        
    def test_haversine_antipodes_large_sphere(self):
        '''Test haversine distance calculation for antipodal points on unit sphere.'''
        first_point = np.array([0,0,self.large_radius])
        second_point = np.array([0,0,-self.large_radius])
        calculated_distance = voronoi_analysis_library.calculate_haversine_distance_between_spherical_points(first_point, second_point, self.large_radius)
        target_distance = self.large_circumference / 2.
        self.assertAlmostEqual(calculated_distance, target_distance, places=7, msg="The haversine distance between antipodal points on the large sphere was not calculated correctly. Target: {target}; calculated: {calc}.".format(target=target_distance, calc=calculated_distance))
        
    def test_haversine_half_arc_unit_sphere(self):
        '''Test haversine distance calculation for half-arc points on unit sphere.'''
        first_point = np.array([0,0,1])
        second_point = np.array([0,1,0])
        calculated_distance = voronoi_analysis_library.calculate_haversine_distance_between_spherical_points(first_point, second_point, self.small_radius)
        target_distance = self.small_circumference / 4.
        self.assertAlmostEqual(calculated_distance, target_distance, places=7, msg="The haversine distance between half-arc points on the unit sphere was not calculated correctly. Target: {target}; calculated: {calc}.".format(target=target_distance, calc=calculated_distance))

    def test_haversine_half_arc_large_sphere(self):
        '''Test haversine distance calculation for half-arc points on large sphere.'''
        first_point = np.array([0,0,self.large_radius])
        second_point = np.array([0,self.large_radius,0])
        calculated_distance = voronoi_analysis_library.calculate_haversine_distance_between_spherical_points(first_point, second_point, self.large_radius)
        target_distance = self.large_circumference / 4.
        self.assertAlmostEqual(calculated_distance, target_distance, places=7, msg="The haversine distance between half-arc points on the large sphere was not calculated correctly. Target: {target}; calculated: {calc}.".format(target=target_distance, calc=calculated_distance))

class TestSphericalCartesianConversion(unittest.TestCase):
    '''Test(s) for conversion of spherical coordinate arrays to cartesian coordinate arrays.'''

    @classmethod
    def setUpClass(cls):
        cls.degrees_input = np.array([[5.0, 20.0, 17.0],
                                      [5.0, 99.3, 66.1]])

        #expected output taken from an online calculator: http://keisan.casio.com/exec/system/1359534351
        cls.expected_degree_output = np.array([[1.373697667, 0.4999850618, 4.78152378],
                                               [-0.7387346631, 4.51118371, 2.025707934]])

        cls.radians_input = np.array([[9.0,1.45,2.66],
                                      [9.0,2.99,0.67]])

        cls.expected_radian_output = np.array([[0.5023424715, 4.138343874, -7.976325095],
                                               [-5.524779676, 0.8439910026, 7.054394993]])


    @classmethod
    def tearDownClass(cls):
        del cls.degrees_input
        del cls.expected_degree_output
        del cls.expected_radian_output
        del cls.radians_input

    def test_degree_input_spherical_to_cartesian(self):
        '''Test conversion of spherical to Cartesian coords using degrees for angles.'''
        actual_result = voronoi_analysis_library.convert_spherical_array_to_cartesian_array(self.degrees_input, angle_measure = 'degrees')
        np.testing.assert_array_almost_equal(actual_result, self.expected_degree_output, decimal = 5, err_msg = "The conversion of spherical coordinates to Cartesian coordinates using degree units for angles does not produce the expected result to the desired precision.")

    def test_radians_input_spherical_to_cartesian(self):
        '''Test conversion of spherical to Cartesian coords using radians for angles.'''
        actual_result = voronoi_analysis_library.convert_spherical_array_to_cartesian_array(self.radians_input, angle_measure = 'radians')
        np.testing.assert_array_almost_equal(actual_result, self.expected_radian_output, decimal = 5, err_msg = "The conversion of spherical coordinates to Cartesian coordinates using radians units for angles does not produce the expected result to the desired precision.")

class TestCartesianSphericalConversion(unittest.TestCase):
    '''Test(s) for conversion of Cartesian coordinate arrays to spherical coordinate arrays.'''

    @classmethod
    def setUpClass(cls):
        cls.Cartesian_input = np.array([[155.6, 178.9, 0.91],
                                        [99.999,3.1015,-8789.2]])

        #expected spherical coord output from an online calculator (http://keisan.casio.com/exec/system/1359533867)
        cls.expected_spherical_output_degrees = np.array([[237.10208371079, 48.984570421492, 89.780097726068],
                                                  [8789.7693973905, 1.7764768911113, 179.34783277354]])

        cls.expected_spherical_output_radians = np.array([[237.10208371079, 0.85494203653007, 1.5669583080822],
                                                         [8789.7693973905, 0.031005370835486, 3.1302101882145]])

    @classmethod
    def tearDownClass(cls):
        del cls.Cartesian_input
        del cls.expected_spherical_output_degrees
        del cls.expected_spherical_output_radians

    def test_cartesian_to_spherical_degrees(self):
        '''Test conversion of a Cartesian coordinate array to a spherical coordinate array with degrees used for angles in the output.'''
        actual_value = voronoi_analysis_library.convert_cartesian_array_to_spherical_array(self.Cartesian_input, angle_measure='degrees')
        np.testing.assert_array_almost_equal(actual_value, self.expected_spherical_output_degrees, decimal = 5, err_msg = "The conversion of cartesian coordinates to spherical coordinates using degree angle measurements did not achieve the expected values to the desired precision.")

    def test_cartesian_to_spherical_radians(self):
        '''Test conversion of a Cartesian coordinate array to a spherical coordinate array with radians used for angles in the output.'''
        actual_value = voronoi_analysis_library.convert_cartesian_array_to_spherical_array(self.Cartesian_input, angle_measure='radians')
        np.testing.assert_array_almost_equal(actual_value, self.expected_spherical_output_radians, decimal = 5, err_msg = "The conversion of cartesian coordinates to spherical coordinates using radian angle measurements did not achieve the expected values to the desired precision.")

class TestRandomSphericalGenerators(unittest.TestCase):
    '''Unit test(s) for the production of random points on the surface of a sphere.'''

    @classmethod
    def setUpClass(cls):
        cls.prng_small = numpy.random.RandomState(1998)
        cls.prng_large = numpy.random.RandomState(556)
        cls.small_num_generators = 55
        cls.large_num_generators = 905
        cls.sphere_radius = 1.1178
        cls.random_generators_small = voronoi_analysis_library.generate_random_array_spherical_generators(cls.small_num_generators, cls.sphere_radius, cls.prng_small)
        cls.random_generators_large = voronoi_analysis_library.generate_random_array_spherical_generators(cls.large_num_generators, cls.sphere_radius, cls.prng_large)

    @classmethod
    def tearDownClass(cls):
        del cls.prng_small
        del cls.prng_large
        del cls.small_num_generators
        del cls.large_num_generators
        del cls.sphere_radius 
        del cls.random_generators_small
        del cls.random_generators_large

    def test_shape_random_spherical_generators(self):
        '''Test the shapes of the data structures produced for small and large random sets of generators on a sphere.'''
        self.assertEqual(self.random_generators_small.shape, (self.small_num_generators, 3), "Incorrect shape for small array of random spherical generators.")
        self.assertEqual(self.random_generators_large.shape, (self.large_num_generators, 3), "Incorrect shape for large array of random spherical generators.")

    def test_radii_random_spherical_systems(self):
        '''Test that the random sets of points on a sphere have the correct radii.'''
        small_distance_matrix = scipy.spatial.distance.pdist(self.random_generators_small)
        large_distance_matrix = scipy.spatial.distance.pdist(self.random_generators_large)
        actual_small_radius = small_distance_matrix.max() / 2.0
        actual_large_radius = large_distance_matrix.max() / 2.0
        self.assertAlmostEqual(actual_small_radius, self.sphere_radius, places = 2, msg="Incorrect radius for small random set of spherical generators. Actual: {actual}; target: {target}".format(actual=actual_small_radius,target=self.sphere_radius))
        self.assertAlmostEqual(actual_large_radius, self.sphere_radius, places = 2, msg="Incorrect radius for large random set of spherical generators. Actual: {actual}; target: {target}".format(actual=actual_small_radius,target=self.sphere_radius))


class TestProduceTrajList(unittest.TestCase):
    '''Unit test(s) for produce_list_trajectories() function.'''

    @classmethod
    def setUpClass(cls):
        cls.d = TempDirectory()
        cls.xtcs = ['dummy1.xtc','dummy2.xtc','dummy3.xtc','dummy4.xtc']
        for filename in cls.xtcs:
            filepath = cls.d.path + '/' + filename
            with open(filepath, 'w') as outfile:
                outfile.write(filename)
        cls.list_of_trajectories = voronoi_analysis_library.produce_list_trajectories(cls.d.path + '/', '*xtc')

    @classmethod
    def tearDownClass(cls):
        cls.d.cleanup()
        del cls.xtcs 
        del cls.list_of_trajectories

    def test_produce_list_trajs(self):
        '''Quick test for proper globbing by produce_list_trajectories().'''
        actual_length = len(self.list_of_trajectories)
        self.assertEqual(actual_length, 4, "produce_list_trajectories is not picking up the full set of xtcs in the dummy directory. Target length: 4; actual length: {length}".format(length=actual_length))
    
    def test_recovered_xtc_files(self):
        '''Test for xtc file names recovered by produce_list_trajectories() function.'''
        xtc_filenames = [filepath.split('/')[-1] for filepath in self.list_of_trajectories]
        self.assertEqual(sorted(xtc_filenames), self.xtcs)

class TestControlCoordProduction(unittest.TestCase):
    '''Unit test(s) for create_control_universe_coord_data() function.'''

    @classmethod
    def setUpClass(cls):
        #need a small flu-like input file (only need residue names to be correct, don't even need coords on sphere for input)
        cls.u = MDAnalysis.Universe('control_traj_2.gro') #using the output of the tested function as the input, so make sure to wipe the coords first
        cls.all_selection = cls.u.select_atoms('all')
        cls.all_selection.set_positions(np.zeros(3)) #all particles are at [0,0,0]
        cls.d = TempDirectory()
        cls.test_input_path = cls.d.path + '/test_input.gro'
        cls.all_selection.atoms.write(cls.test_input_path)
        voronoi_analysis_library.create_control_universe_coord_data(cls.test_input_path, output_path=cls.d.path + '/')
        cls.large_output_universe = MDAnalysis.Universe(cls.d.path + '/' + 'control_traj_1.gro')
        cls.small_output_universe = MDAnalysis.Universe(cls.d.path + '/' + 'control_traj_2.gro')

    @classmethod
    def tearDownClass(cls):
        cls.d.cleanup()
        del cls.all_selection
        del cls.test_input_path
        del cls.large_output_universe 
        del cls.small_output_universe 

    def test_control_coord_radii(self):
        '''Test for the correct radii of inner and outer leaflets in control coordinate data produced.'''
        inner_leaflet_selection_small_system = self.small_output_universe.select_atoms('resname DOPE or resname DOPX or resname POPS')
        inner_leaflet_selection_large_system = self.large_output_universe.select_atoms('resname DOPE or resname DOPX or resname POPS')
        outer_leaflet_selection_small_system = self.small_output_universe.select_atoms('resname CHOL or resname PPCH')
        outer_leaflet_selection_large_system = self.large_output_universe.select_atoms('resname CHOL or resname PPCH')
        inner_small_radius_actual = int(voronoi_analysis_library.convert_cartesian_array_to_spherical_array(inner_leaflet_selection_small_system.coordinates()[0])[0])
        inner_large_radius_actual = int(voronoi_analysis_library.convert_cartesian_array_to_spherical_array(inner_leaflet_selection_large_system.coordinates()[0])[0])
        outer_small_radius_actual = int(voronoi_analysis_library.convert_cartesian_array_to_spherical_array(outer_leaflet_selection_small_system.coordinates()[0])[0])
        outer_large_radius_actual = int(voronoi_analysis_library.convert_cartesian_array_to_spherical_array(outer_leaflet_selection_large_system.coordinates()[0])[0])
        self.assertEqual([inner_small_radius_actual, inner_large_radius_actual, outer_small_radius_actual, outer_large_radius_actual], [500,500,700,700], "The control coordinate leaflet radii are incorrect: {inner_1} {inner_2} {outer_1} {outer_2} instead of 500 500 700 700.".format(inner_1=inner_small_radius_actual,inner_2=inner_large_radius_actual,outer_1=outer_small_radius_actual, outer_2=outer_large_radius_actual))


    def test_num_output_coords(self):
        '''Check that the sizes of the generated control data structures are appropriate.'''
        expected_total_particles = self.all_selection.n_atoms
        upper_limit_total_particles_small = int(0.6 * expected_total_particles) #should be ~half
        actual_particles_small_sys = self.small_output_universe.select_atoms('all').n_atoms
        actual_particles_large_sys = self.large_output_universe.select_atoms('all').n_atoms
        self.assertEqual(actual_particles_large_sys, expected_total_particles, "Total number of particles in large control system is not correct.")
        self.assertLess(actual_particles_small_sys, upper_limit_total_particles_small, "Total number of particles in small control system is not correct (Exceeds 60% of particles in original system).")


class TestControlTrajProduction(unittest.TestCase):
    '''Unit test(s) for create_control_universe_data() function.'''

    @classmethod
    def setUpClass(cls):
        #need a small flu-like input file (only need residue names to be correct, don't even need coords on sphere for input)
        cls.u = MDAnalysis.Universe('control_traj_2.gro') #using the output of the tested function as the input, so make sure to wipe the coords first
        cls.all_selection = cls.u.select_atoms('all')
        cls.all_selection.set_positions(np.zeros(3)) #all particles are at [0,0,0]
        cls.d = TempDirectory()
        cls.test_input_path = cls.d.path + '/test_input.gro'
        cls.all_selection.atoms.write(cls.test_input_path)
        voronoi_analysis_library.create_control_universe_data(cls.test_input_path, output_path=cls.d.path + '/', num_frames=3)
        #call the control *coordinate* version of the function to generate the necessary topology files for reading the trajectories
        voronoi_analysis_library.create_control_universe_coord_data(cls.test_input_path, output_path=cls.d.path + '/') #this function is tested individually elsewhere
        cls.large_output_universe = MDAnalysis.Universe(cls.d.path + '/' + 'control_traj_1.gro', cls.d.path + '/' + 'control_traj_1.xtc')
        cls.small_output_universe = MDAnalysis.Universe(cls.d.path + '/' + 'control_traj_2.gro', cls.d.path + '/' + 'control_traj_2.xtc')

    @classmethod
    def tearDownClass(cls):
        cls.d.cleanup()
        del cls.all_selection
        del cls.test_input_path
        del cls.large_output_universe 
        del cls.small_output_universe 

    def test_control_traj_num_frames(self):
        '''Test that both of the control trajectories produced have the anticipated number of frames.'''
        small_traj_frames = self.small_output_universe.trajectory.n_frames
        large_traj_frames = self.large_output_universe.trajectory.n_frames
        self.assertEqual(small_traj_frames, 3, "Small control trajectory should have 3 frames, but has {actual_frames}.".format(actual_frames=small_traj_frames))
        self.assertEqual(large_traj_frames, 3, "Large control trajectory should have 3 frames, but has {actual_frames}.".format(actual_frames=large_traj_frames))

    def test_similarity_coords_all_frames(self):
        '''Test that the particle coordinates are the same in each frame of both control trajectories (this is a feature of the function).'''
        small_traj = self.small_output_universe.trajectory
        large_traj = self.large_output_universe.trajectory
        small_traj_coords_first_frame = self.small_output_universe.select_atoms('all').coordinates()
        large_traj_coords_first_frame = self.large_output_universe.select_atoms('all').coordinates()
        small_traj[-1]
        large_traj[-1]
        small_traj_coords_final_frame = self.small_output_universe.select_atoms('all').coordinates()
        large_traj_coords_final_frame = self.large_output_universe.select_atoms('all').coordinates()
        np.testing.assert_array_almost_equal(small_traj_coords_first_frame, small_traj_coords_final_frame, decimal = 5, err_msg = "First and final frame coordinates do not match to the desired precision for the small control trajectory.")
        np.testing.assert_array_almost_equal(large_traj_coords_first_frame, large_traj_coords_final_frame, decimal = 5, err_msg = "First and final frame coordinates do not match to the desired precision for the large control trajectory.")


class TestRemoteUniverse(unittest.TestCase):
    '''Unit test(s) for creation of universe objects on remote IPython engines. Probably will not try to create a separate process, but can still make sure universe objects are generated, etc.'''

    @classmethod
    def setUpClass(cls):
        cls.generated_dengue_u = voronoi_analysis_library.produce_universe_object_on_remote_engine_dengue(GRO, XTC)
        #need a temporary directory with some artificial xtc files to properly test the generic universe object generation
        cls.template_u = MDAnalysis.Universe(GRO)
        cls.d = TempDirectory()
        for dummy_xtc_filename in ['dummy_no_solvent_01.xtc','dummy_no_solvent_03.xtc','dummy_no_solvent_05.xtc','dummy_no_solvent_11.xtc', 'dummy_no_solvent_09.xtc']:
            with XTCWriter(cls.d.path + '/' + dummy_xtc_filename, cls.template_u.trajectory.n_atoms) as W:
                W.write(cls.template_u)
        cls.d2 = TempDirectory()
        for dummy_xtc_filename in ['dummy_no_solvent_99.xtc', 'dummy_no_solvent_81.xtc']:
            with XTCWriter(cls.d2.path + '/' + dummy_xtc_filename, cls.template_u.trajectory.n_atoms) as W:
                W.write(cls.template_u)
        cls.generated_u = voronoi_analysis_library.produce_universe_object_on_remote_engine(data_path_1 = cls.d.path + '/', data_path_2 = cls.d2.path + '/', limit_1=-6, limit_2=-4, limit_3=-6,limit_4=-4,coordinate_filepath=GRO, traj_data_extension=1)

    @classmethod
    def tearDownClass(cls):
        cls.d.cleanup()
        cls.d2.cleanup()
        del cls.generated_dengue_u
        del cls.generated_u
        del cls.template_u

    def test_atom_total_generated_dengue_universe(self):
        '''Simply test that artificial dengue universe object has more than 0 atoms.'''
        num_atoms_dengue_u = self.generated_dengue_u.select_atoms('all').n_atoms
        self.assertGreater(num_atoms_dengue_u, 0, "The number of atoms in the artificial dengue universe is not greater than zero.")

    def test_frame_total_generated_dengue_universe(self):
        num_frames_dengue_u = self.generated_dengue_u.trajectory.n_frames
        self.assertGreater(num_frames_dengue_u, 0, "The number of frames in the artificial dengue universe is not greater than zero.")

    def test_frame_total_generated_universe(self):
        num_frames_u = self.generated_u.trajectory.n_frames
        self.assertEqual(num_frames_u, 7, "Incorrect total number of frames in general universe.")

    def test_atom_total_generated_universe(self):
        num_atoms_u = self.generated_u.select_atoms('all').n_atoms
        self.assertGreater(num_atoms_u, 0, "The number of atoms in the artificial general universe is not greater than 0.")

class TestRadialDistanceDengue(unittest.TestCase):
    '''For testing voronoi_analysis_library.precursor_radial_distance_analysis_dengue()'''

    @classmethod
    def setUpClass(cls):
        cls.dengue_universe = MDAnalysis.Universe('dengue_final_snapshot_compact_no_solvent.gro.bz2')
        cls.d = TempDirectory()
        with XTCWriter(cls.d.path + '/' + 'small_dengue_u_testing.xtc', cls.dengue_universe.trajectory.n_atoms) as W:
            for frame in xrange(31): #using 31 frames because of skip 10 in code
                W.write(cls.dengue_universe)
        cls.generated_dengue_universe = MDAnalysis.Universe('dengue_final_snapshot_compact_no_solvent.gro.bz2', cls.d.path + '/small_dengue_u_testing.xtc')
        cls.radial_analysis_output_data_list = voronoi_analysis_library.precursor_radial_distance_analysis_dengue(cls.generated_dengue_universe)

    @classmethod
    def tearDownClass(cls):
        del cls.dengue_universe
        cls.d.cleanup()
        del cls.generated_dengue_universe
        del cls.radial_analysis_output_data_list


    def test_frame_numbers(self):
        '''The list of frame numbers should be [0,10,20,30].'''
        actual_list_frame_numbers = self.radial_analysis_output_data_list[4]
        self.assertEqual(actual_list_frame_numbers, [0,10,20,30], "The list of frame numbers should be [0,10,20,30], but got {actual}".format(actual=actual_list_frame_numbers))
    
    def test_lipid_headgroup_distance_properties(self):
        '''Test some simple properties of the lipid headgroup radial distances.'''
        array_minimum_lipid_headgroup_distances = np.array(self.radial_analysis_output_data_list[0])
        array_maximum_lipid_headgroup_distances = np.array(self.radial_analysis_output_data_list[1])
        array_average_lipid_headgroup_distances = np.array(self.radial_analysis_output_data_list[2])
        array_std_lipid_headgroup_distances = np.array(self.radial_analysis_output_data_list[3])
        array_midpoint_lipid_headgroup_distances = np.array(self.radial_analysis_output_data_list[7])
        np.testing.assert_array_less(array_minimum_lipid_headgroup_distances, array_maximum_lipid_headgroup_distances)
        np.testing.assert_array_less(array_average_lipid_headgroup_distances, array_maximum_lipid_headgroup_distances)
        np.testing.assert_array_less(array_minimum_lipid_headgroup_distances, array_average_lipid_headgroup_distances)
        np.testing.assert_array_less(array_midpoint_lipid_headgroup_distances, array_maximum_lipid_headgroup_distances)
        np.testing.assert_array_less(array_std_lipid_headgroup_distances, array_average_lipid_headgroup_distances)

    def test_headgroup_threshold_percentages(self):
        '''Test to ensure that there are always more headgroups in the outer leaflet than the inner leaflet, for obvious reasons.
        Also, test to ensure that 100 % of headgroups are accounted for.'''
        array_percent_lipid_headgroups_outer_leaflet = np.array(self.radial_analysis_output_data_list[5])
        array_percent_lipid_headgroups_inner_leaflet = np.array(self.radial_analysis_output_data_list[6])
        np.testing.assert_array_less(array_percent_lipid_headgroups_inner_leaflet, array_percent_lipid_headgroups_outer_leaflet)
        np.testing.assert_array_almost_equal(array_percent_lipid_headgroups_inner_leaflet + array_percent_lipid_headgroups_outer_leaflet, numpy.zeros(array_percent_lipid_headgroups_outer_leaflet.shape) + 100., decimal = 1)

    def test_protein_radial_distances(self):
        '''Simple test(s) for the calculated protein radial distances.'''
        array_max_protein_distances = np.array(self.radial_analysis_output_data_list[9])
        array_min_protein_distances = np.array(self.radial_analysis_output_data_list[8])
        np.testing.assert_array_less(array_min_protein_distances, array_max_protein_distances)
        np.testing.assert_array_less(array_max_protein_distances, numpy.zeros(array_max_protein_distances.shape) + 300) #radial distance ceiling check
        np.testing.assert_array_less(numpy.zeros(array_max_protein_distances.shape) + 100, array_max_protein_distances) #radial distance basement check

class TestRadialDistance(unittest.TestCase):
    '''Tests for the general radial distance analysis code -- voronoi_analysis_library.precursor_radial_distance_analysis().'''

    @classmethod
    def setUpClass(cls):
        cls.control_universe = MDAnalysis.Universe('control_traj_2_small.gro.bz2')
        cls.flu_universe_NO_FORS = MDAnalysis.Universe('sim35_final_snapshot_compact_no_solvent_small.gro.bz2')
        cls.flu_universe_WITH_FORS = MDAnalysis.Universe('sim39_final_snapshot_compact_no_solvent_small.gro.bz2')

        cls.d = TempDirectory()
        for temp_testing_xtc_filename, universe in zip(['small_control_testing.xtc', 'small_flu_NO_FORS.xtc','small_flu_WITH_FORS.xtc'],[cls.control_universe, cls.flu_universe_NO_FORS, cls.flu_universe_WITH_FORS]):
            with XTCWriter(cls.d.path + '/' + temp_testing_xtc_filename, universe.trajectory.n_atoms) as W:
                for frame in xrange(301): #using 301 frames because of skip 100 in code
                    W.write(universe)
        cls.generated_control_universe = MDAnalysis.Universe('control_traj_2_small.gro.bz2', cls.d.path + '/small_control_testing.xtc')
        cls.generated_flu_universe_NO_FORS = MDAnalysis.Universe('sim35_final_snapshot_compact_no_solvent_small.gro.bz2', cls.d.path + '/small_flu_NO_FORS.xtc')
        cls.generated_flu_universe_WITH_FORS = MDAnalysis.Universe('sim39_final_snapshot_compact_no_solvent_small.gro.bz2', cls.d.path + '/small_flu_WITH_FORS.xtc')
        cls.radial_analysis_output_data_list_control = voronoi_analysis_library.precursor_radial_distance_analysis(cls.generated_control_universe, FORS_present=None,control_condition=1)
        cls.radial_analysis_output_data_list_flu_NO_FORS = voronoi_analysis_library.precursor_radial_distance_analysis(cls.generated_flu_universe_NO_FORS, FORS_present=None)
        cls.radial_analysis_output_data_list_flu_WITH_FORS = voronoi_analysis_library.precursor_radial_distance_analysis(cls.generated_flu_universe_WITH_FORS,FORS_present=1)

    @classmethod
    def tearDownClass(cls):
        del cls.control_universe
        del cls.flu_universe_NO_FORS
        del cls.flu_universe_WITH_FORS
        cls.d.cleanup()
        del cls.generated_control_universe
        del cls.generated_flu_universe_NO_FORS
        del cls.generated_flu_universe_WITH_FORS
        del cls.radial_analysis_output_data_list_control
        del cls.radial_analysis_output_data_list_flu_NO_FORS
        del cls.radial_analysis_output_data_list_flu_WITH_FORS

    def test_frames_parsed(self):
        '''Ensure that the correct frames are parsed in the radial distance analysis code for each of the three conditions (control, flu, flu + FORS).'''
        actual_list_frame_numbers_control = self.radial_analysis_output_data_list_control[4]
        actual_list_frame_numbers_flu_NO_FORS = self.radial_analysis_output_data_list_flu_NO_FORS[4]
        actual_list_frame_numbers_flu_WITH_FORS = self.radial_analysis_output_data_list_flu_WITH_FORS[4]
        self.assertEqual(actual_list_frame_numbers_control, [0,100,200,300], "The list of control condition frame numbers should be [0,100,200,300], but got {actual}".format(actual=actual_list_frame_numbers_control))
        self.assertEqual(actual_list_frame_numbers_flu_NO_FORS, [0,100,200,300], "The list of flu (NO FORS) condition frame numbers should be [0,100,200,300], but got {actual}".format(actual=actual_list_frame_numbers_flu_NO_FORS))
        self.assertEqual(actual_list_frame_numbers_flu_WITH_FORS, [0,100,200,300], "The list of flu (WITH FORS) condition frame numbers should be [0,100,200,300], but got {actual}".format(actual=actual_list_frame_numbers_flu_WITH_FORS))

    def test_relative_radial_distances(self):
        '''Perform some simple tests of radial distance properties in the three conditions.'''
        PPCH_PO4_min_distances_control = np.array(self.radial_analysis_output_data_list_control[0])
        PPCH_PO4_max_distances_control = np.array(self.radial_analysis_output_data_list_control[1])
        PPCH_PO4_avg_distances_control = np.array(self.radial_analysis_output_data_list_control[2])
        PPCH_PO4_std_distances_control = np.array(self.radial_analysis_output_data_list_control[3])

        PPCH_PO4_min_distances_flu_NO_FORS = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[0])
        PPCH_PO4_max_distances_flu_NO_FORS = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[1])
        PPCH_PO4_avg_distances_flu_NO_FORS = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[2])
        PPCH_PO4_std_distances_flu_NO_FORS = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[3])

        FORS_AM2_min_distances_flu_WITH_FORS = np.array(self.radial_analysis_output_data_list_flu_WITH_FORS[0])
        FORS_AM2_max_distances_flu_WITH_FORS = np.array(self.radial_analysis_output_data_list_flu_WITH_FORS[1])
        FORS_AM2_avg_distances_flu_WITH_FORS = np.array(self.radial_analysis_output_data_list_flu_WITH_FORS[2])
        FORS_AM2_std_distances_flu_WITH_FORS = np.array(self.radial_analysis_output_data_list_flu_WITH_FORS[3])

        #min < max
        np.testing.assert_array_less(PPCH_PO4_min_distances_flu_NO_FORS, PPCH_PO4_max_distances_flu_NO_FORS) 
        np.testing.assert_array_less(FORS_AM2_min_distances_flu_WITH_FORS, FORS_AM2_max_distances_flu_WITH_FORS) 

        #avg < max   
        np.testing.assert_array_less(PPCH_PO4_avg_distances_flu_NO_FORS, PPCH_PO4_max_distances_flu_NO_FORS) 
        np.testing.assert_array_less(FORS_AM2_avg_distances_flu_WITH_FORS, FORS_AM2_max_distances_flu_WITH_FORS) 

        #std < avg
        np.testing.assert_array_less(PPCH_PO4_std_distances_flu_NO_FORS, PPCH_PO4_avg_distances_flu_NO_FORS) 
        np.testing.assert_array_less(FORS_AM2_std_distances_flu_WITH_FORS,FORS_AM2_avg_distances_flu_WITH_FORS) 

        #exact values for control:
        np.testing.assert_array_almost_equal(PPCH_PO4_max_distances_control, numpy.zeros(PPCH_PO4_max_distances_control.shape) + 711.9, decimal = 1)
        np.testing.assert_array_almost_equal(PPCH_PO4_min_distances_control, numpy.zeros(PPCH_PO4_min_distances_control.shape) + 687.9, decimal = 1)

    def test_thresholds(self):
        '''Test some of the leaflet assignment thresholds from radial distance analysis for flu and control conditions.'''
        control_percent_DOPE_DOPX_POPS_below_threshold = np.array(self.radial_analysis_output_data_list_control[19])
        control_percent_DOPE_DOPX_POPS_above_threshold = np.array(self.radial_analysis_output_data_list_control[18])
        control_percent_CHOL_above_threshold = np.array(self.radial_analysis_output_data_list_control[11])

        np.testing.assert_array_almost_equal(control_percent_CHOL_above_threshold, np.zeros(control_percent_CHOL_above_threshold.shape) + 100., decimal=1)
        np.testing.assert_array_less(control_percent_DOPE_DOPX_POPS_above_threshold, control_percent_DOPE_DOPX_POPS_below_threshold) 
        
        flu_NO_FORS_CHOL_above_threshold = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[11])
        flu_NO_FORS_CHOL_below_threshold = np.array(self.radial_analysis_output_data_list_flu_NO_FORS[12])

        np.testing.assert_array_almost_equal(flu_NO_FORS_CHOL_above_threshold + flu_NO_FORS_CHOL_below_threshold, np.zeros(flu_NO_FORS_CHOL_above_threshold.shape) + 100., decimal=1)
        #np.testing.assert_array_less(flu_NO_FORS_CHOL_below_threshold, flu_NO_FORS_CHOL_above_threshold) 

        flu_WITH_FORS_FORS_above_threshold = np.array(self.radial_analysis_output_data_list_flu_WITH_FORS[5])
        np.testing.assert_array_less(np.zeros(flu_WITH_FORS_FORS_above_threshold.shape) + 95.0, flu_WITH_FORS_FORS_above_threshold) #vast majority of FORS should be outer leaflet
