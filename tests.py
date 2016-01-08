import unittest
import numpy as np
import numpy.testing
import voronoi_analysis_library

class TestRawNeighbourAnalysis(unittest.TestCase):

    def setUp(self):
        self.small_data_dict = {'POPC':{},'PPCE':{}} #build a dict with simple coords to test neighbour properties
        self.template_array = np.zeros((3,3))
        self.a = self.template_array[:]
        self.b = self.template_array[:]
        self.c = self.template_array[:]
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

    def test_optimized_neighbour_counting_identical_shuffled_coords(self):
        '''Same test as similarly named method in this class, but for the optimized neighbour search code.'''
        for resname, subdict in self.small_data_dict.iteritems():
            subdict['voronoi_cell_list_vertex_arrays'] = [[self.a,self.b,self.c]]
            subdict['voronoi_cell_list_vertex_arrays_inner_leaflet'] = [[self.a,self.b,self.c]]
        voronoi_neighbour_instance = voronoi_analysis_library.voronoi_neighbour_analysis_optimized(self.small_data_dict, inner_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays_inner_leaflet',outer_leaflet_vertex_list_key = 'voronoi_cell_list_vertex_arrays')
        inner_leaflet_neighbour_dict, outer_leaflet_neighbour_dict = voronoi_neighbour_instance.identify_voronoi_neighbours(frame_index = 0)
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'].keys(), [0], "Inner leaflet POPC should have zero neighbours.") 
        self.assertEqual(inner_leaflet_neighbour_dict['PPCE'].keys(), [0], "Inner leaflet PPCE should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['POPC'].keys(), [0], "Outer leaflet POPC should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'].keys(), [0], "Outer leaflet PPCE should have zero neighbours.") 
        self.assertEqual(outer_leaflet_neighbour_dict['PPCE'][0]['frequency'], 3, "Outer leaflet PPCE should have 3 cells (a,b,c) with 0 neighbours.")
        self.assertEqual(inner_leaflet_neighbour_dict['POPC'][0]['frequency'], 3, "Inner leaflet POPC should have 3 cells (a,b,c) with 0 neighbours.")
