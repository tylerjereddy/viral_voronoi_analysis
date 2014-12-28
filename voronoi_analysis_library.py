'''Library of utility functions for analyzing virus simulations with spherical Voronoi diagrams.'''
import math
import numpy
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle

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

def plot_sample_Voronoi_diagrams(matplotlib_figure_object,list_Voronoi_indices,dict_key_Voronoi_data,plot_title,dict_data):
    plot_number = 1
    color_list = ['black','blue','green','red','purple']
    for current_voronoi_index in list_Voronoi_indices:
        ax = matplotlib_figure_object.add_subplot(1,4,plot_number,projection='3d')
        index = 0
        list_residue_names = []
        for residue_name, subdictionary in dict_data.iteritems():
            list_Voronoi_cell_vertex_arrays = subdictionary[dict_key_Voronoi_data][current_voronoi_index]
            color = color_list[index]
            for vertex_array in list_Voronoi_cell_vertex_arrays:
                polygon = Poly3DCollection([vertex_array/10.],alpha=1.0) #convert to nm
                polygon.set_color(color)
                ax.add_collection3d(polygon)
            list_residue_names.append(residue_name)
            index += 1
        ax.set_title('~{time} $\mu$s ({title})'.format(time=plot_number,title=plot_title))
        ax.auto_scale_xyz
        ax.legend()
        ax.set_xlim(-40,40);ax.set_ylim(-40,40);ax.set_zlim(-40,40);
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')
        ax.set_zlabel('z (nm)')
        POPS_legend_object = Rectangle((0, 0), 1, 1, fc="black")
        DOPE_legend_object = Rectangle((0, 0), 1, 1, fc="blue")
        CHOL_legend_object = Rectangle((0, 0), 1, 1, fc="green")
        PPCH_legend_object = Rectangle((0, 0), 1, 1, fc="red")
        DOPX_legend_object = Rectangle((0, 0), 1, 1, fc="purple")
        ax.legend([POPS_legend_object,DOPE_legend_object,CHOL_legend_object,PPCH_legend_object,DOPX_legend_object],list_residue_names,loc=2,prop={'size':8})
        plot_number += 1
    matplotlib_figure_object.set_size_inches(24,6)


