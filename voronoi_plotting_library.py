'''Library of plotting functions for analysis of virus simulations with spherical Voronoi diagrams.'''

import os
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle
import numpy

class plot_voronoi_neighbour_data_species_specific:
    '''Plot Voronoi neighbour data probing species-specific effects.'''

    def __init__(self,matplotlib_figure_object,inner_leaflet_dict_by_species, outer_leaflet_dict_by_species, simulation_type):
        self.fig = matplotlib_figure_object
        self.subplot_number = 1
        self.inner_leaflet_dict = inner_leaflet_dict_by_species
        self.outer_leaflet_dict = outer_leaflet_dict_by_species
        if simulation_type == 'flu':
            self.color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'} #this should be fine for control, which simply lacks protein
        else:
            self.color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'} 

    def plot(self, num_lipid_species, timestamp_list_microseconds, list_additional_inner_leaflet_dicts, list_additional_outer_leaflet_dicts, area_range, max_num_neighbours = 12, outfile_path = None, hspace = 0.3, wspace = 0.3):
        list_inner_leaflet_dicts = [self.inner_leaflet_dict] + list_additional_inner_leaflet_dicts
        list_outer_leaflet_dicts = [self.outer_leaflet_dict] + list_additional_outer_leaflet_dicts
        current_time_index = 0
        for time in timestamp_list_microseconds:
            inner_leaflet_dict = list_inner_leaflet_dicts[current_time_index]
            outer_leaflet_dict = list_outer_leaflet_dicts[current_time_index]
            for leaflet_name, leaflet_dict in zip(['inner','outer'],[inner_leaflet_dict, outer_leaflet_dict]):
                for lipid_name, neighbours in leaflet_dict.iteritems():
                    ax = self.fig.add_subplot(num_lipid_species * 2, len(timestamp_list_microseconds), self.subplot_number)
                    for neighbour_species_name, neighbour_count_dict in neighbours.iteritems():
                        list_neighbour_counts = []
                        list_avg_surface_areas = []
                        list_std_dev_values = []
                        for neighbour_count, list_surface_areas in neighbour_count_dict.iteritems():
                            list_neighbour_counts.append(neighbour_count)
                            surface_area_array = numpy.array(list_surface_areas)
                            average_surface_area = numpy.average(surface_area_array)
                            std_surface_area = numpy.std(surface_area_array)
                            list_avg_surface_areas.append(average_surface_area)
                            list_std_dev_values.append(std_surface_area)
                        #ax.bar(numpy.array(list_neighbour_counts) - 0.4, list_avg_surface_areas, yerr = numpy.array(list_std_dev_values), alpha = 0.4)
                        ax.errorbar(numpy.array(list_neighbour_counts), list_avg_surface_areas, yerr = None, alpha = 1.0, label = neighbour_species_name, color = self.color_dict[neighbour_species_name])
                        error_array = numpy.array(list_std_dev_values)
                        array_avg_surface_areas = numpy.array(list_avg_surface_areas)
                        ax.fill_between(numpy.array(list_neighbour_counts), array_avg_surface_areas - error_array, array_avg_surface_areas + error_array, alpha = 0.10)

                        ax.set_xlabel('num neighbours', fontsize = 3,labelpad = 0.1)
                        ax.set_ylabel('avg Voronoi cell\n surface area ($\AA^2$)', fontsize = 3, labelpad=0.1)
                        ax.set_xticks(numpy.arange(0,max_num_neighbours, 2))
                        ax.set_title(leaflet_name + ' leaflet ' + lipid_name + ' ({time} $\mu$s)'.format(time = time), fontsize = 3,y = 0.92)
                        ax.tick_params(axis = 'x', labelsize = 3)
                        ax.tick_params(axis = 'y', labelsize = 3)
                        ax.tick_params(which='major', length = 1)
                        lg = ax.legend(prop={'size':1.8}, bbox_to_anchor=(0.97,0.95))
                        fr = lg.get_frame()
                        fr.set_lw(0.2) #thinner legend frame
                        ax.tick_params(direction='in', pad = 1)
                        ax.set_ylim(area_range)
                    self.subplot_number += len(timestamp_list_microseconds)
            current_time_index += 1
            self.subplot_number = current_time_index + 1
        self.fig.set_size_inches(8.5,11) 
        self.fig.subplots_adjust(hspace = hspace, wspace = wspace)
        if outfile_path is not None:
            self.fig.savefig(outfile_path, dpi = 300, bbox_inches = 'tight')

class plot_voronoi_neighbour_data_raw(plot_voronoi_neighbour_data_species_specific):
    '''Plot Voronoi raw neighbour data results.'''

    def plot(self, species_count_dictionary, sim_name, aggregate_figure_object, color_dict, fig_width_inches, fig_height_inches, list_additional_data_dictionaries_inner_leaflet = None, list_additional_data_dictionaries_outer_leaflet = None, list_time_stamps = None, general_plots_xmax = None, aggregate_ymax = None, hspace = 0.4, wspace = 0.5):
        '''species_count_dictionary is for verification that total frequencies of a given molecular species match the amount of that species in the system'''
        full_list_inner_leaflet_data_dicts = [self.inner_leaflet_dict] + list_additional_data_dictionaries_inner_leaflet
        full_list_outer_leaflet_data_dicts = [self.outer_leaflet_dict] + list_additional_data_dictionaries_outer_leaflet
        assert len(full_list_inner_leaflet_data_dicts) == len(full_list_outer_leaflet_data_dicts), "Inner and Outer leaflet data dictionary lists must have the same length."
        num_columns = 2 * len(full_list_inner_leaflet_data_dicts)
        num_rows = 2 * len(species_count_dictionary.keys())
        num_rows_aggregate = 2
        aggregate_subplot_counter = 1
        list_subplot_numbers_current_time_stamp = []
        for row in xrange(num_rows):
            list_subplot_numbers_current_time_stamp.append(numpy.array([1,2]) + (row * num_columns))



        #print list_subplot_numbers_current_time_stamp
        for inner_leaflet_dict, outer_leaflet_dict, time_stamp in zip(full_list_inner_leaflet_data_dicts, full_list_outer_leaflet_data_dicts, list_time_stamps): #iterate over multiple simulation time points
            ax_aggregate_inner_area = aggregate_figure_object.add_subplot(2,num_columns,aggregate_subplot_counter)
            ax_aggregate_inner_freq = aggregate_figure_object.add_subplot(2,num_columns,aggregate_subplot_counter + 1)
            ax_aggregate_outer_area = aggregate_figure_object.add_subplot(2,num_columns,aggregate_subplot_counter + num_columns)
            ax_aggregate_outer_freq = aggregate_figure_object.add_subplot(2,num_columns,aggregate_subplot_counter + num_columns + 1)
            list_aggregate_frequency_data_inner_leaflet = [] #if 5 neighbours occurs 70 times overall, then want 70 copies of number 5 in this list, etc. (to facilitate histogram data structure)
            list_aggregate_frequency_data_outer_leaflet = []
            subplot_index = 0
            for molecular_species_name in species_count_dictionary.keys():
                inner_species_counter = 0
                outer_species_counter = 0
                frequency_count_molecular_species = 0
                for leaflet_name, leaflet_dict in zip(['inner','outer'],[inner_leaflet_dict, outer_leaflet_dict]):
                    #print molecular_species_name, leaflet_name, time_stamp, list_subplot_numbers_current_time_stamp[subplot_index]
                    #print molecular_species_name, leaflet_name, time_stamp, list_subplot_numbers_current_time_stamp[subplot_index]
                    ax = self.fig.add_subplot(num_rows,num_columns,list_subplot_numbers_current_time_stamp[subplot_index][0])
                    ax2 = self.fig.add_subplot(num_rows,num_columns,list_subplot_numbers_current_time_stamp[subplot_index][1])
                    for axis in [ax, ax2]:
                        axis.tick_params(axis = 'x', labelsize = 3)
                        axis.tick_params(axis = 'y', labelsize = 3)
                        axis.tick_params(which='major', length = 1)
                    num_neighbours_list = []
                    frequency_of_neighbour_count_list = []
                    for num_neighbours, subdictionary in leaflet_dict[molecular_species_name].iteritems():
                        frequency_of_neighbour_count = subdictionary['frequency']
                        list_surface_area_for_voronoi_cells_with_this_neighbour_count = subdictionary['list_surface_areas']
                        array_surface_areas = numpy.array(list_surface_area_for_voronoi_cells_with_this_neighbour_count)
                        average_surface_area = numpy.average(array_surface_areas)
                        std_surface_area = numpy.std(array_surface_areas)
                        ax.errorbar(num_neighbours,average_surface_area, yerr = std_surface_area, marker = 'o',c='blue', markeredgecolor='none', ms = 2, capsize=1, elinewidth=0.5)
                        ax.set_xlabel('# neighbours', fontsize =3, labelpad = 0.1)
                        ax.set_ylabel('avg Voronoi cell\n surface area ($\AA^2$)', fontsize =3, labelpad = 0.1)
                        if general_plots_xmax:
                            ax.set_xlim(-0.2,general_plots_xmax)
                        else:
                            ax.set_xlim(-0.2,12)
                            
                        num_neighbours_list.append(num_neighbours)
                        frequency_of_neighbour_count_list.append(frequency_of_neighbour_count)
                        if leaflet_name == 'inner':
                            list_aggregate_frequency_data_inner_leaflet.extend(frequency_of_neighbour_count * [num_neighbours])
                            if inner_species_counter == 0:
                                ax_aggregate_inner_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name],label=molecular_species_name)
                                ax_aggregate_inner_area.set_title('{condition} aggregate inner Leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp), fontsize =3)
                                inner_species_counter += 1
                            else:
                                ax_aggregate_inner_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name])
                        else:
                            list_aggregate_frequency_data_outer_leaflet.extend(frequency_of_neighbour_count * [num_neighbours])
                            if outer_species_counter == 0:
                                ax_aggregate_outer_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name],label=molecular_species_name)
                                ax_aggregate_outer_area.set_title('{condition} aggregate outer Leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp), fontsize =3)
                                outer_species_counter += 1
                            else:
                                ax_aggregate_outer_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name])
                                
                    ax2.bar(numpy.array(num_neighbours_list) - 0.3, frequency_of_neighbour_count_list, lw = 0)
                    if general_plots_xmax:
                        ax2.set_xlim(-1,general_plots_xmax)
                    else:
                        ax2.set_xlim(-1,12)
                    ax2.set_ylabel('frequency', fontsize =3, labelpad = 0.1)
                    ax2.set_xlabel('# neighbours', fontsize =3, labelpad = 0.1)
                    for axis in [ax,ax2]:
                        axis.set_title('{condition} {leaflet_flag}\nleaflet {species} ({time_stamp})'.format(leaflet_flag = leaflet_name, species = molecular_species_name, condition = sim_name, time_stamp = time_stamp), fontsize = 3, y = 0.92)
                        axis.tick_params(direction='in', pad = 1)
                    subplot_index += 1
                    #print 'subplot_index:', subplot_index
                    frequency_count_molecular_species += sum(frequency_of_neighbour_count_list)
                assert frequency_count_molecular_species == species_count_dictionary[molecular_species_name], "The neighbour frequency count for {mol_species} does not match the total molecules of this type in the system. Got {actual_count} instead.".format(mol_species = molecular_species_name, actual_count = frequency_count_molecular_species)
                
            for axis in [ax_aggregate_inner_area, ax_aggregate_outer_area]:
                axis.legend(loc=2,fontsize=3,scatterpoints=1)
                axis.set_xlabel('# neighbours', fontsize = 3)
                axis.set_ylabel('avg Voronoi cell\n surface area ($\AA^2$)', fontsize = 3, labelpad = 0.1)
                if general_plots_xmax:
                    axis.set_xlim(0,general_plots_xmax)
                else:
                    axis.set_xlim(0,14)
                axis.set_ylim(0,350)
                
            for axis in [ax_aggregate_inner_freq,ax_aggregate_outer_freq]:
                axis.set_xlabel('number of Voronoi neighbours', fontsize = 3)
                if general_plots_xmax:
                    axis.set_xlim(0,general_plots_xmax)
                else:
                    axis.set_xlim(0,14)
                axis.set_ylabel('Frequency', fontsize = 3)
                if aggregate_ymax:
                    axis.set_ylim(0,aggregate_ymax)
                else:
                    axis.set_ylim(0,2000)
                
            assert len(list_aggregate_frequency_data_inner_leaflet) + len(list_aggregate_frequency_data_outer_leaflet) == sum(species_count_dictionary.itervalues()), "The total number of Voronoi cells for which there are requency values should match the total number of particles in the system."

            bins = numpy.arange(14) + 0.1
            ax_aggregate_inner_freq.set_title('{condition} aggregate inner leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp), fontsize =3)
            ax_aggregate_inner_freq.hist(list_aggregate_frequency_data_inner_leaflet,bins = bins)
            ax_aggregate_outer_freq.set_title('{condition} aggregate outer leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp), fontsize =3)
            ax_aggregate_outer_freq.hist(list_aggregate_frequency_data_outer_leaflet,bins = bins)


            list_subplot_numbers_current_time_stamp = numpy.array(list_subplot_numbers_current_time_stamp) + 2
            aggregate_subplot_counter += 2








#            ax = self.fig.add_subplot(100,2,self.subplot_number)
#            ax.set_xlabel('num neighbours')
#            ax.set_ylabel('avg Voronoi cell\n surface area ($\AA^2$)')
#            ax.set_title(leaflet_name + ' leaflet ' + lipid_name)
#            self.subplot_number += 1
        self.fig.set_size_inches(fig_width_inches,fig_height_inches) 
        self.fig.subplots_adjust(hspace = hspace, wspace = wspace)
        aggregate_figure_object.set_size_inches(8,3) 
        #aggregate_figure_object(hspace = 0.7, wspace = 0.3)

def plot_sample_Voronoi_diagrams(matplotlib_figure_object,list_Voronoi_indices,dict_key_Voronoi_data,plot_title,dict_data,dengue_condition=None,control_condition=None):
    plot_number = 1

    if not dengue_condition:
        color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'} #this should be fine for control, which simply lacks protein
    else:
        color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'} 

    for current_voronoi_index in list_Voronoi_indices:
        ax = matplotlib_figure_object.add_subplot(1,4,plot_number,projection='3d')
        index = 0
        list_residue_names = []
        for residue_name, subdictionary in dict_data.iteritems():
            list_Voronoi_cell_vertex_arrays = subdictionary[dict_key_Voronoi_data][current_voronoi_index]
            color = color_dict[residue_name]
            for vertex_array in list_Voronoi_cell_vertex_arrays: #control condition debug -- indexing 0:10 produces some tiny polygons, 0:100 causes a straight kernel crash
                #ax.plot(vertex_array[...,0]/10.,vertex_array[...,1]/10.,vertex_array[...,2]/10.,c=color,linestyle='',marker='.',markeredgecolor='none') #debug scatter for control condition
                polygon = Poly3DCollection([vertex_array/10.],alpha=1.0) #convert to nm
                polygon.set_color(color)
                ax.add_collection3d(polygon)
            list_residue_names.append(residue_name)
            index += 1
        ax.set_title('~{time} $\mu$s ({title})'.format(time=plot_number,title=plot_title))
        ax.auto_scale_xyz
        ax.legend()
        if not dengue_condition:
            if not control_condition:
                ax.set_xlim(-40,40);ax.set_ylim(-40,40);ax.set_zlim(-40,40);
            else: #the control condition
                ax.set_xlim(-80,80);ax.set_ylim(-80,80);ax.set_zlim(-80,80);
        else:
            ax.set_xlim(-25,25);ax.set_ylim(-25,25);ax.set_zlim(-25,25);
        ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')
        ax.set_zlabel('z (nm)')
        list_legend_objects = [Rectangle((0, 0), 1, 1, fc=color_dict[residue_name]) for residue_name in list_residue_names]
        ax.legend(list_legend_objects,list_residue_names,loc=2,prop={'size':8})
        plot_number += 1
    matplotlib_figure_object.set_size_inches(24,6)

def plot_sample_Voronoi_diagrams_zoom(matplotlib_figure_object,list_Voronoi_indices,dict_key_Voronoi_data,plot_title,dict_data,dengue_condition=None,debug_generators=None,generator_array=None, outfile_name=None):
    plot_number = 1
    if not dengue_condition:
        color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'}
    else:
        color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'} 
    if debug_generators:
        alpha_value = 0.5
    else:
        alpha_value = 1.0
    for current_voronoi_index in list_Voronoi_indices:
        ax = matplotlib_figure_object.add_subplot(1,4,plot_number,projection='3d')
        index = 0
        list_residue_names = []
        for residue_name, subdictionary in dict_data.iteritems():
            list_Voronoi_cell_vertex_arrays = subdictionary[dict_key_Voronoi_data][current_voronoi_index]
            color = color_dict[residue_name]
            for vertex_array in list_Voronoi_cell_vertex_arrays:
                if numpy.abs(vertex_array[...,1:]).max() > 99: #filter beyond zoom limit
                    continue
                if vertex_array[...,0].max() < 0: #filter beyond zoom limit
                    continue
                polygon = Poly3DCollection([vertex_array/10.],alpha=alpha_value) #convert to nm
                polygon.set_color(color)
                polygon.set_edgecolor('black')
                ax.add_collection3d(polygon)
            list_residue_names.append(residue_name)
            index += 1
        ax.set_title('~{time} $\mu$s ({title})'.format(time=plot_number,title=plot_title))
        #it is proving tricky to place labels and ticks in this 3D zoom-in workflow
        #ax.text(-10,0,0,'x (nm)',fontsize=60)
        #ax.set_xlabel('x (nm)')
        ax.set_ylabel('y (nm)')
        ax.set_zlabel('z (nm)')
        list_legend_objects = [Rectangle((0, 0), 1, 1, fc=color_dict[residue_name]) for residue_name in list_residue_names]
        #ax.legend(list_legend_objects,list_residue_names,loc=2,prop={'size':8})
        ax.set_xlim(0,40);
        ax.set_ylim(-10,10);ax.set_zlim(-10,10);
        ax.w_xaxis.set_ticklabels([''])
        ax.azim = 0
        ax.elev = 0
        if debug_generators:
            if 'outer' in plot_title:
                filtered_generator_array = generator_array[(generator_array[...,0] > 650) & (numpy.abs(generator_array[...,1]) <= 99) & (numpy.abs(generator_array[...,2]) <= 99)]
            else:
                filtered_generator_array = generator_array[(generator_array[...,0] < 600) & (generator_array[...,0] > 0) & (numpy.abs(generator_array[...,1]) <= 99) & (numpy.abs(generator_array[...,2]) <= 99)]
            ax.plot(filtered_generator_array[...,0]/10.,filtered_generator_array[...,1]/10.,filtered_generator_array[...,2]/10.,c='k',marker='.',linestyle='')
        plot_number += 1
    matplotlib_figure_object.set_size_inches(24,6)
    #matplotlib_figure_object.savefig('control_1_inner_leaflet_pydata_london_2015.png', dpi = 300)
    if outfile_name is not None:
        matplotlib_figure_object.savefig(outfile_name, dpi = 300, bbox_inches = 'tight')

def area_per_molecule_plotting(figure_object,list_frame_numbers,list_percent_surface_area_reconstitution=None,list_percent_surface_area_reconstitution_inner_leaflet=None,protein_present=None,simulation_title=None,dictionary_headgroup_data=None,list_percent_surface_area_reconstitution_from_lipids_only=None,list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet=None,list_percent_surface_area_reconstitution_from_proteins_only=None,list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet=None,control_condition=None,dengue_condition=None, output_file=None):
    if not dengue_condition:
        color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'}
    else: #dengue has a very different lipid profile
        color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'}
    if not protein_present:
        ax_inner = figure_object.add_subplot('141')
        ax_outer = figure_object.add_subplot('142')
        if not control_condition:
            array_time_values = numpy.array(list_frame_numbers) / 10000. #microseconds
        else:
            array_time_values = numpy.array(list_frame_numbers) / 1000. #microseconds
        array_percent_surface_area_reconstitution = numpy.array(list_percent_surface_area_reconstitution)
        ax_outer.scatter(array_time_values,array_percent_surface_area_reconstitution,c='black',edgecolor='None',label='outer leaflet')
        array_percent_surface_area_reconstitution_inner_leaflet = numpy.array(list_percent_surface_area_reconstitution_inner_leaflet)
        ax_inner.scatter(array_time_values,array_percent_surface_area_reconstitution_inner_leaflet,c='red',edgecolor='None',label='inner leaflet')
        for axis in [ax_outer,ax_inner]:
            axis.set_ylim(90,101)
            axis.set_xlim(0,5)
            axis.legend(loc=4)
            axis.set_ylabel('Percent surface area reconstitution\n from Voronoi cells')
            axis.set_xlabel('Time ($\mu$s)')
            axis.set_title(simulation_title)

        ax2 = figure_object.add_subplot('143')
        index = 0
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            color = color_dict[residue_name]
            array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list'])
            array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list'])
            ax2.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
            #ax2.fill_between(array_frame_numbers,array_voronoi_cell_areas-array_voronoi_cell_std_dev,array_voronoi_cell_areas+array_voronoi_cell_std_dev,color=color,alpha=0.2) 
            index += 1
#ax2.legend(loc=4)
        ax2.set_ylabel('Average area per lipid ($\AA^2$)')
        ax2.set_xlim(0,5)
        ax2.set_xlabel('Time ($\mu$s)')
        ax2.set_title(simulation_title + ' outer leaflet')
        if not control_condition:
            ax2.set_ylim(34,54)
        else:
            ax2.set_ylim(180,550)
            pass

        ax4 = figure_object.add_subplot('144')
        index = 0
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            color = color_dict[residue_name]
            array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list_inner_leaflet'])
            array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list_inner_leaflet'])
            ax4.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
            #ax2.fill_between(array_frame_numbers,array_voronoi_cell_areas-array_voronoi_cell_std_dev,array_voronoi_cell_areas+array_voronoi_cell_std_dev,color=color,alpha=0.2) 
            index += 1
        ax4.legend()
        ax4.set_ylabel('Average area per lipid ($\AA^2$)')
        ax4.set_xlim(0,5)
        ax4.set_xlabel('Time ($\mu$s)')
        ax4.set_title(simulation_title + ' inner leaflet')
        if not control_condition:
            ax4.set_ylim(34,54)
        else:
            ax4.set_ylim(180,550)


    else: #protein is present
        ax_1_inner = figure_object.add_subplot('141')
        ax_1_outer = figure_object.add_subplot('142')
        if not dengue_condition:
            array_time_values = numpy.array(list_frame_numbers) / 10000. #microseconds
        else:
            array_time_values = numpy.array(list_frame_numbers) / 1000. #microseconds
        array_percent_surface_area_reconstitution_lipid_outer_leaflet = numpy.array(list_percent_surface_area_reconstitution_from_lipids_only)
        array_percent_surface_area_reconstitution_lipid_inner_leaflet = numpy.array(list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet)
        array_percent_surface_area_reconstitution_protein_outer_leaflet = numpy.array(list_percent_surface_area_reconstitution_from_proteins_only)
        array_percent_surface_area_reconstitution_protein_inner_leaflet = numpy.array(list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet)
        #print out values for citation in dengue paper
        print '---------'
        print simulation_title
        print 'avg array_percent_surface_area_reconstitution_protein_inner_leaflet:', numpy.average(array_percent_surface_area_reconstitution_protein_inner_leaflet)
        print 'avg array_percent_surface_area_reconstitution_protein_outer_leaflet:', numpy.average(array_percent_surface_area_reconstitution_protein_outer_leaflet)
        print '---------'
        combined_percent_reconstitution_array_outer_leaflet = array_percent_surface_area_reconstitution_lipid_outer_leaflet + array_percent_surface_area_reconstitution_protein_outer_leaflet
        combined_percent_reconstitution_array_inner_leaflet = array_percent_surface_area_reconstitution_lipid_inner_leaflet + array_percent_surface_area_reconstitution_protein_inner_leaflet

        ax_1_outer.scatter(array_time_values,array_percent_surface_area_reconstitution_lipid_outer_leaflet,c='black',edgecolor='None',label='lipid outer',marker='o')
        ax_1_outer.scatter(array_time_values,array_percent_surface_area_reconstitution_protein_outer_leaflet,c='black',edgecolor='None',label='protein outer',marker='^')
        ax_1_outer.scatter(array_time_values,combined_percent_reconstitution_array_outer_leaflet,c='black',edgecolor='None',label='combined outer',marker='*',s=50)
        ax_1_inner.scatter(array_time_values,array_percent_surface_area_reconstitution_lipid_inner_leaflet,c='red',edgecolor='None',label='lipid inner',marker='o',alpha=0.5)
        ax_1_inner.scatter(array_time_values,array_percent_surface_area_reconstitution_protein_inner_leaflet,c='red',edgecolor='None',label='protein inner',marker='^',alpha=0.5)
        ax_1_inner.scatter(array_time_values,combined_percent_reconstitution_array_inner_leaflet,c='red',edgecolor='None',label='combined inner',marker='*',s=50,alpha=0.5)
        for axis in [ax_1_outer, ax_1_inner]:
            axis.set_ylim(-10,110)
            axis.set_xlim(0,5)
            axis.set_ylabel('Percent surface area reconstitution\n from Voronoi cells')
            axis.set_xlabel('Time ($\mu$s)')
            axis.legend(loc=0,bbox_to_anchor=[1.0, 0.5],ncol=2,fontsize=8)
            axis.set_title(simulation_title)

        ax_2 = figure_object.add_subplot('143')
        index = 0
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            color = color_dict[residue_name]
            array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list'])
            array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list'])
            ax_2.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
            index += 1
#ax_2.legend(loc=1)
        ax_2.set_ylabel('Average area per molecule ($\AA^2$)')
        ax_2.set_xlim(0,5)
        ax_2.set_xlabel('Time ($\mu$s)')
        ax_2.set_title(simulation_title + ' outer leaflet')
        ax_2.set_ylim(20,200)
#ax_2.set_ylim(34,54)

        ax_4 = figure_object.add_subplot('144')
        index = 0
        for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
            color = color_dict[residue_name]
            array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list_inner_leaflet'])
            array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list_inner_leaflet'])
            ax_4.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
            index += 1
        ax_4.legend(loc=0)
        ax_4.set_ylabel('Average area per molecule ($\AA^2$)')
        ax_4.set_xlim(0,5)
        ax_4.set_ylim(20,200)
        ax_4.set_xlabel('Time ($\mu$s)')
        ax_4.set_title(simulation_title + ' inner leaflet')
#ax_4.set_ylim(34,54)

    if output_file is not None:
        figure_object.subplots_adjust(wspace = 0.5)
        figure_object.set_size_inches(8.5,3)

        for axis in [ax_outer, ax_inner, ax2, ax4]:
            axis.set_title(axis.get_title(), fontsize=6)
            if axis == ax_outer or axis == ax_inner:
                axis.legend(loc = 4, prop={'size':6}, scatterpoints=1)
            else:
                if 'control_2' in output_file:
                    location = 4
                else:
                    location = 0
                axis.legend(loc = location, prop={'size':6}, scatterpoints=1)
            ylbl = axis.yaxis.get_label()
            old_ylbl_fontsize = ylbl.get_fontsize()
            ylbl.set_fontsize(6)
                    

        figure_object.savefig(output_file, dpi = 300, bbox_inches = 'tight')

        #reset values after producing images so I don't throw off the render in the Jupyter nb
        for axis in [ax_outer, ax_inner, ax2, ax4]:
            axis.set_title(axis.get_title(), fontdict=None)
            if axis == ax_outer or axis == ax_inner:
                axis.legend(loc = 4, prop={'size':10})
            else:
                axis.legend(loc = 0, prop={'size':10})
            ylbl = axis.yaxis.get_label()
            ylbl.set_fontsize(old_ylbl_fontsize)

    figure_object.subplots_adjust(wspace = 0.2)
    figure_object.set_size_inches(25,4)

def area_per_molecule_plotting_sysL_figures(figure_object_1,figure_object_2,figure_name_1,figure_name_2,list_frame_numbers,list_percent_surface_area_reconstitution=None,list_percent_surface_area_reconstitution_inner_leaflet=None,protein_present=None,simulation_title=None,dictionary_headgroup_data=None,list_percent_surface_area_reconstitution_from_lipids_only=None,list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet=None,list_percent_surface_area_reconstitution_from_proteins_only=None,list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet=None,control_condition=None,dengue_condition=None):
    if not dengue_condition:
        color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'}
        array_time_values = numpy.array(list_frame_numbers) / 10000. #microseconds
    else: #dengue has a very different lipid profile
        color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'}
        array_time_values = numpy.array(list_frame_numbers) / 1000. #microseconds

    ax2 = figure_object_1.add_subplot('111')
    index = 0
    for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
        color = color_dict[residue_name]
        array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list'])
        array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list'])
        ax2.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
        #ax2.fill_between(array_frame_numbers,array_voronoi_cell_areas-array_voronoi_cell_std_dev,array_voronoi_cell_areas+array_voronoi_cell_std_dev,color=color,alpha=0.2) 
        index += 1
#ax2.legend(loc=4)
    ax2.set_ylabel('Avg area / molecule ($\AA^2$)')
    ax2.set_xlim(0,5)
    ax2.set_xlabel('Time ($\mu$s)')
    #ax2.set_title(simulation_title + ' outer leaflet')
    ax2.set_ylim(20,160)

    ax4 = figure_object_2.add_subplot('111')
    index = 0
    for residue_name, subdictionary in dictionary_headgroup_data.iteritems():
        color = color_dict[residue_name]
        array_voronoi_cell_areas = numpy.array(subdictionary['voronoi_cell_avg_values_list_inner_leaflet'])
        array_voronoi_cell_std_dev = numpy.array(subdictionary['voronoi_cell_std_values_list_inner_leaflet'])
        ax4.scatter(array_time_values,array_voronoi_cell_areas,label=residue_name,edgecolor='None',color=color)
        #ax2.fill_between(array_frame_numbers,array_voronoi_cell_areas-array_voronoi_cell_std_dev,array_voronoi_cell_areas+array_voronoi_cell_std_dev,color=color,alpha=0.2) 
        index += 1
    #ax4.legend()
    ax4.set_ylabel('Avg area / molecule ($\AA^2$)')
    ax4.set_xlim(0,5)
    ax4.set_xlabel('Time ($\mu$s)')
    #ax4.set_title(simulation_title + ' inner leaflet')
    ax4.set_ylim(20,160)

    figure_object_1.set_size_inches(2,2)
    figure_object_2.set_size_inches(2,2)
    for fig in [figure_object_1,figure_object_2]:
        fig.subplots_adjust(left=0.3,bottom=0.22)

    figure_object_1.savefig(figure_name_1,dpi=300)
    figure_object_2.savefig(figure_name_2,dpi=300)

class radial_distance_assessment_dengue:

    def __init__(self,matplotlib_figure_object,list_frame_numbers_dengue,list_min_dengue_lipid_headgroup_distances,list_max_dengue_lipid_headgroup_distances,list_average_dengue_lipid_headgroup_distances,list_std_dev_dengue_lipid_headgroup_distances,list_dengue_lipid_headgroup_midpoint_distances,list_dengue_lipid_headgroup_percent_above_threshold,list_dengue_lipid_headgroup_percent_below_threshold,list_min_protein_distances,list_max_protein_distances, array_min_dengue_lipid_headgroup_radial_distances,array_max_dengue_lipid_headgroup_radial_distances,array_min_protein_distances,array_max_protein_distances,array_average_dengue_lipid_headgroup_radial_distances,array_std_dev_dengue_lipid_headgroup_radial_distances,array_dengue_lipid_headgroup_unbiased_midpoint_distances,array_dengue_lipid_headgroup_percent_above_midpoint_threshold,array_dengue_lipid_headgroup_percent_below_midpoint_threshold,array_frame_numbers):

        self.matplotlib_figure_object = matplotlib_figure_object
        #dengue data structure initialization:
        self.array_min_dengue_lipid_headgroup_radial_distances = array_min_dengue_lipid_headgroup_radial_distances
        self.array_max_dengue_lipid_headgroup_radial_distances = array_max_dengue_lipid_headgroup_radial_distances
        self.array_min_protein_distances = array_min_protein_distances
        self.array_max_protein_distances = array_max_protein_distances
        self.array_average_dengue_lipid_headgroup_radial_distances = array_average_dengue_lipid_headgroup_radial_distances
        self.array_std_dev_dengue_lipid_headgroup_radial_distances = array_std_dev_dengue_lipid_headgroup_radial_distances
        self.array_dengue_lipid_headgroup_unbiased_midpoint_distances = array_dengue_lipid_headgroup_unbiased_midpoint_distances
        self.array_dengue_lipid_headgroup_percent_above_midpoint_threshold = array_dengue_lipid_headgroup_percent_above_midpoint_threshold
        self.array_dengue_lipid_headgroup_percent_below_midpoint_threshold = array_dengue_lipid_headgroup_percent_below_midpoint_threshold
        self.array_frame_numbers = array_frame_numbers

    def plot(self,title_string):
        '''Plot the dengue radial distance assessment data.'''
        ax = self.matplotlib_figure_object.add_subplot('121')
        #print self.array_frame_numbers.shape, self.array_min_dengue_lipid_headgroup_radial_distances.shape #debug
        ax.scatter(self.array_frame_numbers,self.array_min_dengue_lipid_headgroup_radial_distances,label='min dengue lipid headgroup radial distance',c='black',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_min_protein_distances,label='min protein distances',c='grey',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_max_protein_distances,label='max protein distances',c='green',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_max_dengue_lipid_headgroup_radial_distances,label='max dengue lipid headgroup radial distance',c='red',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_average_dengue_lipid_headgroup_radial_distances,label='average dengue lipid headgroup radial distance',c='blue',edgecolor='None')
        ax.fill_between(self.array_frame_numbers,self.array_average_dengue_lipid_headgroup_radial_distances-self.array_std_dev_dengue_lipid_headgroup_radial_distances,self.array_average_dengue_lipid_headgroup_radial_distances+self.array_std_dev_dengue_lipid_headgroup_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean dengue lipid headgroup OD values
        ax.scatter(self.array_frame_numbers,self.array_dengue_lipid_headgroup_unbiased_midpoint_distances,label='unbiased dengue lipid headgroup radial midpoints',c='yellow',edgecolor='None')
        ax.set_xlim(-100,5100)
        ax.set_ylim(10,70)
        ax.set_xlabel('Frame #')
        ax.set_ylabel('Radial distance from virion centroid (nm)')
        ax.legend(loc=2,fontsize=8)

        ax2 = self.matplotlib_figure_object.add_subplot('122')
        ax2.scatter(self.array_frame_numbers,self.array_dengue_lipid_headgroup_percent_above_midpoint_threshold,label='above midpoint',color='orange')
        ax2.scatter(self.array_frame_numbers,self.array_dengue_lipid_headgroup_percent_below_midpoint_threshold,label='below midpoint',color='blue')
        ax2.set_ylabel('Percent dengue lipid headgroup particles above\n or below midpoint')
        ax2.set_xlabel('Frame #')
        ax2.set_xlim(-100,5100)
        ax2.legend(loc=2)

        for axis in [ax,ax2]:
            axis.set_title(title_string)

        self.matplotlib_figure_object.set_size_inches(16,6)

class radial_distance_assessment:

    def __init__(self,matplotlib_figure_object,list_min_PPCH_PO4_distances,list_max_PPCH_PO4_distances,list_average_PPCH_PO4_distances,list_std_dev_PPCH_PO4_distances,list_frame_numbers,list_PPCH_percent_above_threshold,list_min_CHOL_ROH_distances,list_max_CHOL_ROH_distances,list_average_CHOL_ROH_distances,list_std_dev_CHOL_ROH_distances,list_CHOL_ROH_midpoint_distances,list_CHOL_ROH_percent_above_threshold,list_CHOL_ROH_percent_below_threshold,list_min_remaining_headgroup_distances,list_max_remaining_headgroup_distances,list_average_remaining_headgroup_distances,list_std_dev_remaining_headgroup_distances,list_remaining_headgroup_midpoint_distances,list_remaining_headgroup_percent_above_threshold,list_remaining_headgroup_percent_below_threshold,PPCH_threshold,list_min_FORS_AM2_distances=None,list_max_FORS_AM2_distances=None,list_average_FORS_AM2_distances=None,list_std_dev_FORS_AM2_distances=None,list_FORS_percent_above_treshold=None,FORS_present=None,control_condition=None):

        self.matplotlib_figure_object = matplotlib_figure_object
        self.threshold = PPCH_threshold
        self.FORS_present = FORS_present
        #PPCH data initialization:
        self.array_min_PPCH_PO4_radial_distances = numpy.array(list_min_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_max_PPCH_PO4_radial_distances = numpy.array(list_max_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_average_PPCH_PO4_radial_distances = numpy.array(list_average_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_std_dev_PPCH_PO4_radial_distances = numpy.array(list_std_dev_PPCH_PO4_distances) / 10.0 #convert to nm
        self.array_percent_PPCH_PO4_above_threshold = numpy.array(list_PPCH_percent_above_threshold)
        if self.FORS_present:
            #FORS data initialization:
            self.array_min_FORS_AM2_radial_distances = numpy.array(list_min_FORS_AM2_distances) / 10.0 #convert to nm
            self.array_max_FORS_AM2_radial_distances = numpy.array(list_max_FORS_AM2_distances) / 10.0 #convert to nm
            self.array_average_FORS_AM2_radial_distances = numpy.array(list_average_FORS_AM2_distances) / 10.0 #convert to nm
            self.array_std_dev_FORS_AM2_radial_distances = numpy.array(list_std_dev_FORS_AM2_distances) / 10.0 #convert to nm
            self.array_percent_FORS_AM2_above_threshold = numpy.array(list_FORS_percent_above_treshold)
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
        self.control_condition = control_condition

    def plot(self,title_string,equil_line=None):
        '''Plot the radial distance assessment data.'''
        ax = self.matplotlib_figure_object.add_subplot('421')
        ax.scatter(self.array_frame_numbers,self.array_min_PPCH_PO4_radial_distances,label='min PPCH PO4 radial distance',c='black',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_max_PPCH_PO4_radial_distances,label='max PPCH PO4 radial distance',c='red',edgecolor='None')
        ax.scatter(self.array_frame_numbers,self.array_average_PPCH_PO4_radial_distances,label='average PPCH PO4 radial distance',c='blue',edgecolor='None')
        ax.fill_between(self.array_frame_numbers,self.array_average_PPCH_PO4_radial_distances-self.array_std_dev_PPCH_PO4_radial_distances,self.array_average_PPCH_PO4_radial_distances+self.array_std_dev_PPCH_PO4_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean PPCH PO4 OD values
        ax.set_xlabel('Frame #')
        ax.set_ylabel('Radial distance from vesicle centroid (nm)')
        if equil_line:
            ax.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
#now, use a second plot to track the % of PPCH PO4 particles that fall above the assigned radial distance threshold
        if not self.control_condition:
            ax.axhline(y=self.threshold/10.,xmin=0,xmax=50000,c='green') #radial distance values above this threshold should capture most of the PPCH PO4 particles (within 1 std dev of the mean)
            ax.set_ylim(20,45)
            ax.set_xlim(-900,50000)
            ax.legend()
        else:
            ax.axhline(y=60,xmin=0,xmax=50000,c='green') #use 60 nm as the cutoff for control conditions, which use absolute leaflet assignments (1 species = 1 leaflet only)
            ax.set_ylim(40,80)
            ax.set_xlim(-90,5000)
            ax.legend(fontsize=8)
        ax2 = self.matplotlib_figure_object.add_subplot('422')
        #print 'self.array_frame_numbers.shape:', self.array_frame_numbers.shape, 'self.array_percent_PPCH_PO4_above_threshold.shape:', self.array_percent_PPCH_PO4_above_threshold.shape #debug
        ax2.scatter(self.array_frame_numbers,self.array_percent_PPCH_PO4_above_threshold,color='orange',edgecolor='None')
        ax2.set_xlabel('Frame #')
        ax2.set_ylabel('Percent PPCH PO4 particles above cutoff\n radial distance threshold')
        ax2.axhline(y=98.0,xmin=0,xmax=50000,c='purple',lw=6,alpha=0.4) #98% of PPCH PO4 particles
        if equil_line:
            ax2.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        if not self.control_condition:
            ax2.set_ylim(80,100.0)
            ax2.set_xlim(-900,50000)
        else:
            ax2.set_ylim(80,105.0)
            ax2.set_xlim(-90,5000)


#now, CHOL-related plots in the second row

        ax3 = self.matplotlib_figure_object.add_subplot('423')
        ax3.scatter(self.array_frame_numbers,self.array_min_CHOL_ROH_radial_distances,label='min CHOL ROH radial distance',c='black',edgecolor='None')
        ax3.scatter(self.array_frame_numbers,self.array_max_CHOL_ROH_radial_distances,label='max CHOL ROH radial distance',c='red',edgecolor='None')
        ax3.scatter(self.array_frame_numbers,self.array_average_CHOL_ROH_radial_distances,label='average CHOL ROH radial distance',c='blue',edgecolor='None')
        ax3.fill_between(self.array_frame_numbers,self.array_average_CHOL_ROH_radial_distances-self.array_std_dev_CHOL_ROH_radial_distances,self.array_average_CHOL_ROH_radial_distances+self.array_std_dev_CHOL_ROH_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean CHOL ROH OD values
        ax3.scatter(self.array_frame_numbers,self.array_CHOL_ROH_unbiased_midpoint_distances,label='unbiased CHOL ROH radial midpoints',c='yellow',edgecolor='None')
        if equil_line:
            ax3.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        if not self.control_condition:
            ax3.set_xlim(-900,50000)
            ax3.set_ylim(20,45)
            ax3.legend()
        else: 
            ax3.set_xlim(-90,5000)
            ax3.set_ylim(40,80)
            ax3.axhline(y=60,xmin=0,xmax=50000,c='green') #use 60 nm as the cutoff for control conditions, which use absolute leaflet assignments (1 species = 1 leaflet only)
            ax3.legend(fontsize=8)
        ax3.set_xlabel('Frame #')
        ax3.set_ylabel('Radial distance from vesicle centroid (nm)')
        ax4 = self.matplotlib_figure_object.add_subplot('424')
        if equil_line:
            ax4.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax4.set_xlabel('Frame #')
        if not self.control_condition:
            ax4.set_xlim(-900,50000)
            ax4.scatter(self.array_frame_numbers,self.array_CHOL_ROH_percent_above_midpoint_threshold,label='above midpoint',color='orange')
            ax4.scatter(self.array_frame_numbers,self.array_CHOL_ROH_percent_below_midpoint_threshold,label='below midpoint',color='blue')
            ax4.set_ylabel('Percent CHOL ROH particles above\n or below midpoint')
            ax4.legend()
        else:
            ax4.set_ylabel('Percent CHOL ROH particles above\nradial distance threshold')
            ax4.scatter(self.array_frame_numbers,self.array_CHOL_ROH_percent_above_midpoint_threshold,color='orange')
            ax4.set_xlim(-90,5000)
            ax4.set_ylim(-5,105)


        ax5 = self.matplotlib_figure_object.add_subplot('425')
        ax5.scatter(self.array_frame_numbers,self.array_min_remaining_headgroup_radial_distances,label='min [DOPE/X, POPS] PO4 radial distance',c='black',edgecolor='None')
        ax5.scatter(self.array_frame_numbers,self.array_max_remaining_headgroup_radial_distances,label='max [DOPE/X, POPS] PO4 radial distance',c='red',edgecolor='None')
        ax5.scatter(self.array_frame_numbers,self.array_average_remaining_headgroup_radial_distances,label='average [DOPE/X, POPS] PO4 radial distance',c='blue',edgecolor='None')
        ax5.fill_between(self.array_frame_numbers,self.array_average_remaining_headgroup_radial_distances-self.array_std_dev_remaining_headgroup_radial_distances,self.array_average_remaining_headgroup_radial_distances+self.array_std_dev_remaining_headgroup_radial_distances,color='blue',alpha=0.2) 
        ax5.scatter(self.array_frame_numbers,self.array_remaining_headgroup_unbiased_midpoint_distances,label='unbiased [DOPE/X, POPS] PO4 radial midpoints',c='yellow',edgecolor='None')
        ax5.set_ylabel('Radial distance from vesicle centroid (nm)')
        if not self.control_condition:
            ax5.set_xlim(-900,50000)
            ax5.set_ylim(20,45)
            ax5.legend()
        else: 
            ax5.set_xlim(-90,5000)
            ax5.set_ylim(40,80)
            ax5.axhline(y=60,xmin=0,xmax=50000,c='green') #use 60 nm as the cutoff for control conditions, which use absolute leaflet assignments (1 species = 1 leaflet only)
            ax5.legend(fontsize=8)
            #ax5.set_ylim(20,45)
        ax5.set_xlabel('Frame #')
        if equil_line:
            ax5.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable
        ax6 = self.matplotlib_figure_object.add_subplot('426')
        ax6.set_xlabel('Frame #')
        if not self.control_condition:
            ax6.set_ylabel('Percent [DOPE/X, POPS] PO4 particles above\n or below midpoint')
            ax6.set_xlim(-900,50000)
            ax6.scatter(self.array_frame_numbers,self.array_remaining_headgroup_percent_above_midpoint_threshold,label='above midpoint',color='orange')
            ax6.scatter(self.array_frame_numbers,self.array_remaining_headgroup_percent_below_midpoint_threshold,label='below midpoint',color='blue')
            ax6.legend()
        else:
            ax6.set_ylabel('Percent [DOPE/X, POPS] PO4 particles above\n radial distance threshold')
            ax6.scatter(self.array_frame_numbers,self.array_remaining_headgroup_percent_above_midpoint_threshold,color='orange')
            ax6.set_xlim(-90,5000)
            ax6.set_ylim(-5,105)
        if equil_line: #300 ns equil line
            ax6.axvline(x=3000,ymin=0,ymax=1,c='green') #300 ns (3000 frame) equilibration line -- where the holes have sealed and the OD is stable

        #plot FORS data, if applicable:
        if self.FORS_present:
            ax7 = self.matplotlib_figure_object.add_subplot('427')
            ax7.scatter(self.array_frame_numbers,self.array_min_FORS_AM2_radial_distances,label='min FORS AM2 radial distance',c='black',edgecolor='None')
            ax7.scatter(self.array_frame_numbers,self.array_max_FORS_AM2_radial_distances,label='max FORS AM2 radial distance',c='red',edgecolor='None')
            ax7.scatter(self.array_frame_numbers,self.array_average_FORS_AM2_radial_distances,label='average FORS AM2 radial distance',c='blue',edgecolor='None')
            ax7.fill_between(self.array_frame_numbers,self.array_average_FORS_AM2_radial_distances-self.array_std_dev_FORS_AM2_radial_distances,self.array_average_FORS_AM2_radial_distances+self.array_std_dev_FORS_AM2_radial_distances,color='blue',alpha=0.2) #show the standard deviation about the mean FORS AM2 OD values
            ax7.set_xlabel('Frame #')
            ax7.set_ylabel('Radial distance from vesicle centroid (nm)')
            ax7.legend()
            ax7.axhline(y=self.threshold/10.,xmin=0,xmax=50000,c='green') #radial distance values above this threshold should capture most of the FORS AM2 particles (within 1 std dev of the mean)
#now, use a second plot to track the % of FORS AM2 particles that fall above the assigned radial distance threshold
            ax7.set_ylim(20,45)
            ax7.set_xlim(-900,50000)
            ax8 = self.matplotlib_figure_object.add_subplot('428')
            #print 'self.array_frame_numbers.shape:', self.array_frame_numbers.shape, 'self.array_percent_FORS_AM2_above_threshold.shape:', self.array_percent_FORS_AM2_above_threshold.shape #debug
            ax8.scatter(self.array_frame_numbers,self.array_percent_FORS_AM2_above_threshold,color='orange',edgecolor='None')
            ax8.set_xlabel('Frame #')
            ax8.set_ylabel('Percent FORS AM2 particles above cutoff\n radial distance threshold')
            ax8.axhline(y=98.0,xmin=0,xmax=50000,c='purple',lw=6,alpha=0.4) #98% of FORS AM2 particles
            ax8.set_ylim(80,100.0)
            ax8.set_xlim(-900,50000)

        if self.FORS_present:
            axis_list = [ax,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
        else:
            axis_list = [ax,ax2,ax3,ax4,ax5,ax6]

        for axis in axis_list:
            axis.set_title(title_string)

        self.matplotlib_figure_object.set_size_inches(16,24)

class plot_sample_N_neighbours(object):
    '''Produce sample zoom-in Voronoi diagrams with a central molecular species surrounded by N neighbours (and nothing else).'''

    def __init__(self, sample_neighbour_dict, name_central_species, condition=None):
        self.central_species_name = name_central_species
        self.central_Voronoi_cell_array = sample_neighbour_dict['central_cell']
        self.neighbour_dict = sample_neighbour_dict['neighbours']
        self.figure = matplotlib.pyplot.figure()
        self.condition = condition
        if condition == 'flu':
            self.azimuth = 180
            self.color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'}
        elif condition == 'dengue':
            self.azimuth = 0
            self.color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'}
        self.list_residue_names = [self.central_species_name]

    def plot(self, plot_title):
        #infer plot limits from Voronoi cell data
        list_x_value_arrays_neighbours = self.central_Voronoi_cell_array[...,0].tolist()
        list_y_value_arrays_neighbours = self.central_Voronoi_cell_array[...,1].tolist()
        list_z_value_arrays_neighbours = self.central_Voronoi_cell_array[...,2].tolist()

        for species_name, list_arrays_Voronoi_cells in self.neighbour_dict.iteritems():
            for vertex_array in list_arrays_Voronoi_cells: 
                list_x_value_arrays_neighbours.extend(vertex_array[...,0])
                list_y_value_arrays_neighbours.extend(vertex_array[...,1])
                list_z_value_arrays_neighbours.extend(vertex_array[...,2])

        x_values = numpy.array(list_x_value_arrays_neighbours)
        y_values = numpy.array(list_y_value_arrays_neighbours)
        z_values = numpy.array(list_z_value_arrays_neighbours)
        
        #equalize axis limits for visualization purposes:
        x_size = abs(x_values.max() - x_values.min())
        y_size = abs(y_values.max() - y_values.min())
        z_size = abs(z_values.max() - z_values.min())
        max_width = numpy.array([x_size, y_size, z_size]).max()
        x_delta = max_width - x_size
        y_delta = max_width - y_size
        z_delta = max_width - z_size
        #balance the widths by adding half the data to the min and max positions:
        x_min = x_values.min() - x_delta / 2.0
        y_min = y_values.min() - y_delta / 2.0
        z_min = z_values.min() - z_delta / 2.0
        x_max = x_values.max() + x_delta / 2.0
        y_max = y_values.max() + y_delta / 2.0
        z_max = z_values.max() + z_delta / 2.0


        #prepare for plotting
        ax = self.figure.add_subplot('111', projection = '3d')
        #ax.set_title(plot_title)

        #plot the central Voronoi cell
        polygon = Poly3DCollection([self.central_Voronoi_cell_array/10.],alpha=0.5) #convert to nm
        color = self.color_dict[self.central_species_name]
        polygon.set_color(color)
        ax.add_collection3d(polygon)
        
        #plot the neighbouring Voronoi cells   
        for species_name, list_arrays_Voronoi_cells in self.neighbour_dict.iteritems():
            if species_name not in self.list_residue_names:
                self.list_residue_names.append(species_name)
            for vertex_array in list_arrays_Voronoi_cells: 
                polygon = Poly3DCollection([vertex_array/10.],alpha=0.5) #convert to nm
                color = self.color_dict[species_name]
                polygon.set_color(color)
                ax.add_collection3d(polygon)

        ax.set_xlim(x_min/10., x_max/10.)
        ax.set_ylim(y_min/10., y_max/10.)
        ax.set_zlim(z_min/10., z_max/10.)
        #ax.set_xlabel('x (nm)', fontsize=16)
        ax.set_ylabel('\ny (nm)', fontsize=16)
        ax.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
        ax.azim = self.azimuth
        ax.elev = 0
        ax.tick_params(axis='both', which='major', labelsize=11)
        list_legend_objects = [Rectangle((0, 0), 1, 1, fc=self.color_dict[residue_name], alpha=0.5) for residue_name in self.list_residue_names]
        if self.condition == 'flu':
            ax.legend(list_legend_objects,self.list_residue_names,loc='center',prop={'size':10}, bbox_to_anchor=[0.2, 0.5])
            ax.set_zlabel('\nz (nm)', fontsize=16)
        elif self.condition == 'dengue':
            ax.legend(list_legend_objects,self.list_residue_names,loc='center',prop={'size':10}, bbox_to_anchor=[0.85, 0.7])
            ax.set_zlabel('z (nm)\n', fontsize=16)
        self.figure.set_size_inches(5,5)

    def save_plot(self, filename):
        self.figure.savefig(filename, dpi = 300)

def plot_panel_fig_isolated_raw(dictionary_data_list, lipid_name,  plot_title, plot_filename, time_value_list_microseconds):
    '''Plot one of the panels for figure isolated_raw.'''
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot('111')
    legend_added = []
    for dictionary_data, colour, time_value in zip(dictionary_data_list, ['red', 'green', 'grey'], time_value_list_microseconds):
        for neighbour_count, subdict in dictionary_data[lipid_name].iteritems():
            avg_surface_area = numpy.average(subdict['list_surface_areas'])
            std_surface_area = numpy.std(subdict['list_surface_areas'])
            if time_value in legend_added:
                ax.errorbar(neighbour_count,avg_surface_area, yerr = std_surface_area, marker = 'o', c = colour, markeredgecolor='none')
            else:
                ax.errorbar(neighbour_count,avg_surface_area, yerr = std_surface_area, marker = 'o', c = colour, markeredgecolor='none', label = str(time_value) + ' $\mu$s')
                legend_added.append(time_value)

    try:
        os.mkdir('./fig_isolated_raw')
    except:
        pass

    ax.set_title(plot_title)
    ax.set_xlim(0,15)
    ax.set_ylim(0,300)
    ax.set_xlabel('Raw Neighbour Count')
    ax.set_ylabel('Avg Surface Area ($\AA^2$)')
    ax.legend(loc=2, numpoints=1)
    fig.set_size_inches(4.3,3.5)
    fig.subplots_adjust(left=0.2)
    fig.savefig(os.path.join('./fig_isolated_raw',plot_filename), dpi = 300)

def plot_panel_fig_aggregate_raw(dictionary_data_list, plot_title, plot_filename, time_value_list_microseconds):
    '''Plot one of the panels for figure aggregate_raw.'''
    fig = matplotlib.pyplot.figure()
    subplots = ['131','132','133']
    if 'Flu' in plot_title:
        color_dict = {'POPS':'black','DOPE':'blue','CHOL':'green','PPCH':'red','DOPX':'purple','protein':'orange','FORS':'brown'}
    else:
        color_dict = {'POPC':'black','PPCE':'blue','DPPE':'green','CER':'red','DUPC':'purple','protein':'orange','DOPS':'brown','PPCS':'pink'}
    for dictionary_data, time_value in zip(dictionary_data_list, time_value_list_microseconds):
        legend_added = []
        ax = fig.add_subplot(subplots.pop(0))
        ax.set_ylim(0,300)
        ax.set_xlim(0,15)
        ax.set_title(plot_title + ' ' + str(time_value) + ' $\mu$s')
        ax.set_xlabel('Raw Neighbour Count')
        ax.set_ylabel('Surface Area ($\AA^2$)')
        ax.set_yticks([0,100,200,300])
        for lipid_name in dictionary_data.keys():
            for neighbour_count, subdict in dictionary_data[lipid_name].iteritems():
                avg_surface_area = numpy.average(subdict['list_surface_areas'])
                std_surface_area = numpy.std(subdict['list_surface_areas'])
                if lipid_name in legend_added:
                    ax.scatter(neighbour_count,avg_surface_area, marker = 'o', c = color_dict[lipid_name], edgecolor='none')
                else:
                    ax.scatter(neighbour_count,avg_surface_area, marker = 'o', c = color_dict[lipid_name], edgecolor='none', label = lipid_name)
                    legend_added.append(lipid_name)

    ax.legend(loc=4, scatterpoints=1, bbox_to_anchor = [1.9, -0.3], ncol=2)

    try:
        os.mkdir('./fig_aggregate_raw')
    except:
        pass

    fig.set_size_inches(13.5, 1.5)
    fig.subplots_adjust(left=0.05, wspace = 0.25, top = 0.72, bottom = 0.29, right = 0.8)
    fig.savefig(os.path.join('./fig_aggregate_raw',plot_filename), dpi = 300)

