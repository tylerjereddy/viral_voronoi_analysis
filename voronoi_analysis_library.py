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

class radial_distance_assessment:

    def __init__(self,matplotlib_figure_object,list_min_PPCH_PO4_distances,list_max_PPCH_PO4_distances,list_average_PPCH_PO4_distances,list_std_dev_PPCH_PO4_distances,list_frame_numbers,list_PPCH_percent_above_threshold,list_min_CHOL_ROH_distances,list_max_CHOL_ROH_distances,list_average_CHOL_ROH_distances,list_std_dev_CHOL_ROH_distances,list_CHOL_ROH_midpoint_distances,list_CHOL_ROH_percent_above_threshold,list_CHOL_ROH_percent_below_threshold,list_min_remaining_headgroup_distances,list_max_remaining_headgroup_distances,list_average_remaining_headgroup_distances,list_std_dev_remaining_headgroup_distances,list_remaining_headgroup_midpoint_distances,list_remaining_headgroup_percent_above_threshold,list_remaining_headgroup_percent_below_threshold,PPCH_threshold):

        self.matplotlib_figure_object = matplotlib_figure_object
        self.threshold = PPCH_threshold
        #PPCH data initialization:
        self.array_min_PPCH_PO4_radial_distances = numpy.array(list_min_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_max_PPCH_PO4_radial_distances = numpy.array(list_max_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_average_PPCH_PO4_radial_distances = numpy.array(list_average_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_std_dev_PPCH_PO4_radial_distances = numpy.array(list_std_dev_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_percent_PPCH_PO4_above_threshold = numpy.array(list_PPCH_percent_above_threshold)
        #CHOL data initialization:
        self.array_min_CHOL_ROH_radial_distances = numpy.array(list_min_CHOL_ROH_distances)/ 10.
        self.array_max_CHOL_ROH_radial_distances = numpy.array(list_max_CHOL_ROH_distances)/ 10.
        self.array_average_CHOL_ROH_radial_distances = numpy.array(list_average_CHOL_ROH_distances)/ 10.
        self.array_std_dev_CHOL_ROH_radial_distances = numpy.array(list_std_dev_CHOL_ROH_distances)/ 10.
        self.array_CHOL_ROH_unbiased_midpoint_distances = numpy.array(list_CHOL_ROH_midpoint_distances) / 10.
        self.array_CHOL_ROH_percent_above_midpoint_threshold = numpy.array(list_CHOL_ROH_percent_above_threshold)
        self.array_CHOL_ROH_percent_below_midpoint_threshold = numpy.array(list_CHOL_ROH_percent_below_threshold)
        #now for the remaining headgroup particles (POPS, DOPX/E)
        self.array_min_remaining_headgroup_radial_distances = numpy.array(list_min_remaining_headgroup_distances) / 10.
        self.array_max_remaining_headgroup_radial_distances = numpy.array(list_max_remaining_headgroup_distances) / 10.
        self.array_average_remaining_headgroup_radial_distances = numpy.array(list_average_remaining_headgroup_distances) / 10.
        self.array_std_dev_remaining_headgroup_radial_distances = numpy.array(list_std_dev_remaining_headgroup_distances) / 10.
        self.array_remaining_headgroup_unbiased_midpoint_distances = numpy.array(list_remaining_headgroup_midpoint_distances) / 10.
        self.array_remaining_headgroup_percent_above_midpoint_threshold = numpy.array(list_remaining_headgroup_percent_above_threshold)
        self.array_remaining_headgroup_percent_below_midpoint_threshold = numpy.array(list_remaining_headgroup_percent_below_threshold)

        self.array_frame_numbers = numpy.array(list_frame_numbers)

    def plot(self,title_string,equil_line=None):
        '''Plot the radial distance assessment data.'''
        ax = self.matplotlib_figure_object.add_subplot('321')
        ax.scatter(self.array_frame_numbers,self.array_min_PPCH_PO4_radial_distances,label='min PPCH PO4 radial distance',c='black',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_max_PPCH_PO4_radial_distances,label='max PPCH PO4 radial distance',c='red',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_average_PPCH_PO4_radial_distances,label='average PPCH PO4 radial distance',c='blue',edgecolor='None')
        ax.fill_between(self.array_frame_numbers,self.array_average_PPCH_PO4_radial_distances-self.array_std_dev_PPCH_PO4_radial_distances,self.array_average_PPCH_PO4_radial_distances+self.array_std_dev_PPCH_PO4_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean PPCH PO4 OD values
        ax.set_xlabel('Frame #')
        ax.set_ylabel('Radial distance from vesicle centroid (nm)')
        ax.legend()
        ax.axhline(y=self.threshold/10.,xmin=0,xmax=50000,c='green') #radial distance values above this threshold should capture most of the PPCH PO4 particles (within 1 std dev of the mean)
        if equil_line:
            ax.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
#now, use a second plot to track the % of PPCH PO4 particles that fall above the assigned radial distance threshold
        ax.set_ylim(20,45)
        ax.set_xlim(-900,50000)
        ax2 = self.matplotlib_figure_object.add_subplot('322')
        ax2.scatter(self.array_frame_numbers,self.array_percent_PPCH_PO4_above_threshold,color='orange',edgecolor='None')
        ax2.set_xlabel('Frame #')
        ax2.set_ylabel('Percent PPCH PO4 particles above cutoff\n radial distance threshold')
        ax2.axhline(y=98.0,xmin=0,xmax=50000,c='purple',lw=6,alpha=0.4) #98% of PPCH PO4 particles
        if equil_line:
            ax2.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax2.set_ylim(80,100.0)
        ax2.set_xlim(-900,50000)

#now, CHOL-related plots in the second row

        ax3 = self.matplotlib_figure_object.add_subplot('323')
        ax3.scatter(self.array_frame_numbers,self.array_min_CHOL_ROH_radial_distances,label='min CHOL ROH radial distance',c='black',edgecolor='None')
        ax3.scatter(self.array_frame_numbers,self.array_max_CHOL_ROH_radial_distances,label='max CHOL ROH radial distance',c='red',edgecolor='None')
        ax3.scatter(self.array_frame_numbers,self.array_average_CHOL_ROH_radial_distances,label='average CHOL ROH radial distance',c='blue',edgecolor='None')
        ax3.fill_between(self.array_frame_numbers,self.array_average_CHOL_ROH_radial_distances-self.array_std_dev_CHOL_ROH_radial_distances,self.array_average_CHOL_ROH_radial_distances+self.array_std_dev_CHOL_ROH_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean CHOL ROH OD values
        ax3.scatter(self.array_frame_numbers,self.array_CHOL_ROH_unbiased_midpoint_distances,label='unbiased CHOL ROH radial midpoints',c='yellow',edgecolor='None')
        if equil_line:
            ax3.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax3.set_xlim(-900,50000)
        ax3.set_ylim(20,45)
        ax3.set_xlabel('Frame #')
        ax3.set_ylabel('Radial distance from vesicle centroid (nm)')
        ax3.legend()
        ax4 = self.matplotlib_figure_object.add_subplot('324')
        ax4.scatter(self.array_frame_numbers,self.array_CHOL_ROH_percent_above_midpoint_threshold,label='above midpoint',color='orange')
        ax4.scatter(self.array_frame_numbers,self.array_CHOL_ROH_percent_below_midpoint_threshold,label='below midpoint',color='blue')
        if equil_line:
            ax4.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax4.set_ylabel('Percent CHOL ROH particles above\n or below midpoint')
        ax4.set_xlabel('Frame #')
        ax4.set_xlim(-900,50000)
        ax4.legend()


        ax5 = self.matplotlib_figure_object.add_subplot('325')
        ax5.scatter(self.array_frame_numbers,self.array_min_remaining_headgroup_radial_distances,label='min [DOPE/X, POPS] PO4 radial distance',c='black',edgecolor='None')
        ax5.scatter(self.array_frame_numbers,self.array_max_remaining_headgroup_radial_distances,label='max [DOPE/X, POPS] PO4 radial distance',c='red',edgecolor='None')
        ax5.scatter(self.array_frame_numbers,self.array_average_remaining_headgroup_radial_distances,label='average [DOPE/X, POPS] PO4 radial distance',c='blue',edgecolor='None')
        ax5.fill_between(self.array_frame_numbers,self.array_average_remaining_headgroup_radial_distances-self.array_std_dev_remaining_headgroup_radial_distances,self.array_average_remaining_headgroup_radial_distances+self.array_std_dev_remaining_headgroup_radial_distances,color='blue',alpha=0.2) 
        ax5.scatter(self.array_frame_numbers,self.array_remaining_headgroup_unbiased_midpoint_distances,label='unbiased [DOPE/X, POPS] PO4 radial midpoints',c='yellow',edgecolor='None')
        ax5.set_ylabel('Radial distance from vesicle centroid (nm)')
        ax5.set_xlim(-900,50000)
        ax5.set_xlabel('Frame #')
        if equil_line:
            ax5.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax5.legend()
        ax5.set_ylim(20,45)
        ax6 = self.matplotlib_figure_object.add_subplot('326')
        ax6.scatter(self.array_frame_numbers,self.array_remaining_headgroup_percent_above_midpoint_threshold,label='above midpoint',color='orange')
        ax6.scatter(self.array_frame_numbers,self.array_remaining_headgroup_percent_below_midpoint_threshold,label='below midpoint',color='blue')
        ax6.set_ylabel('Percent [DOPE/X, POPS] PO4 particles above\n or below midpoint')
        ax6.set_xlabel('Frame #')
        ax6.set_xlim(-900,50000)
        if equil_line: #300 ns equil line
            ax6.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax6.legend()

        for axis in [ax,ax2,ax3,ax4,ax5,ax6]:
            axis.set_title(title_string)

        self.matplotlib_figure_object.set_size_inches(16,24)

def TMD_particle_selector(input_array,molecule_type):
    '''selects the TMD coordinate elements from the input array and combines to a simplified new numpy array with TMD particle coordinates only.'''
    if molecule_type == 'lipid':
        output_array = input_array #still using the COG of the entire lipid
    elif molecule_type == 'HA': #the index numbers are based on study of topology combined with Danny's DPhil thesis
        HA_1_TMD_numpy_array = input_array[1110:1174] 
        HA_2_TMD_numpy_array = input_array[2305:2369]
        HA_3_TMD_numpy_array = input_array[3500:3564]
        output_array = numpy.concatenate((HA_1_TMD_numpy_array,HA_2_TMD_numpy_array,HA_3_TMD_numpy_array))
    elif molecule_type == 'NA':
        NA_1_TMD_numpy_array = input_array[13:70]
        NA_2_TMD_numpy_array = input_array[1037:1094]
        NA_3_TMD_numpy_array = input_array[2059:2116]
        NA_4_TMD_numpy_array = input_array[3083:3140]
        output_array = numpy.concatenate((NA_1_TMD_numpy_array,NA_2_TMD_numpy_array,NA_3_TMD_numpy_array,NA_4_TMD_numpy_array))
    elif molecule_type == 'M2':
        output_array = input_array #this protein doesn't really have an ectodomain so I think it is quite acceptable to continue with usage of the overall COG
    return output_array

