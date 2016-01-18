import unittest
import numpy as np
import numpy.testing
import voronoi_analysis_library

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
