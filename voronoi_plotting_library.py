'''Library of plotting functions for analysis of virus simulations with spherical Voronoi diagrams.'''

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle

class plot_voronoi_neighbour_data_species_specific:
    '''Plot Voronoi neighbour data probing species-specific effects.'''

    def __init__(self,matplotlib_figure_object,inner_leaflet_dict_by_species, outer_leaflet_dict_by_species):
        self.fig = matplotlib_figure_object
        self.subplot_number = 1
        self.inner_leaflet_dict = inner_leaflet_dict_by_species
        self.outer_leaflet_dict = outer_leaflet_dict_by_species

    def plot(self, num_lipid_species, timestamp_list_microseconds, list_additional_inner_leaflet_dicts, list_additional_outer_leaflet_dicts, area_range, max_num_neighbours = 12):
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
                        ax.errorbar(numpy.array(list_neighbour_counts), list_avg_surface_areas, yerr = None, alpha = 1.0, label = neighbour_species_name)
                        error_array = numpy.array(list_std_dev_values)
                        array_avg_surface_areas = numpy.array(list_avg_surface_areas)
                        ax.fill_between(numpy.array(list_neighbour_counts), array_avg_surface_areas - error_array, array_avg_surface_areas + error_array, alpha = 0.05)

                        ax.set_xlabel('num neighbours')
                        ax.set_ylabel('avg Voronoi cell surface area ($\AA^2$)')
                        ax.set_xticks(numpy.arange(0,max_num_neighbours, 2))
                        ax.set_title(leaflet_name + ' leaflet ' + lipid_name + ' ({time} $\mu$s)'.format(time = time))
                        ax.legend(prop={'size':8})
                        ax.set_ylim(area_range)
                    self.subplot_number += len(timestamp_list_microseconds)
            current_time_index += 1
            self.subplot_number = current_time_index + 1
        self.fig.set_size_inches(25,80) 
        self.fig.subplots_adjust(hspace = 0.3, wspace = 0.3)

class plot_voronoi_neighbour_data_raw(plot_voronoi_neighbour_data_species_specific):
    '''Plot Voronoi raw neighbour data results.'''

    def plot(self, species_count_dictionary, sim_name, aggregate_figure_object, color_dict, fig_width_inches, fig_height_inches, list_additional_data_dictionaries_inner_leaflet = None, list_additional_data_dictionaries_outer_leaflet = None, list_time_stamps = None, general_plots_xmax = None, aggregate_ymax = None):
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
                    num_neighbours_list = []
                    frequency_of_neighbour_count_list = []
                    for num_neighbours, subdictionary in leaflet_dict[molecular_species_name].iteritems():
                        frequency_of_neighbour_count = subdictionary['frequency']
                        list_surface_area_for_voronoi_cells_with_this_neighbour_count = subdictionary['list_surface_areas']
                        array_surface_areas = numpy.array(list_surface_area_for_voronoi_cells_with_this_neighbour_count)
                        average_surface_area = numpy.average(array_surface_areas)
                        std_surface_area = numpy.std(array_surface_areas)
                        ax.errorbar(num_neighbours,average_surface_area, yerr = std_surface_area, marker = 'o',c='blue')
                        ax.set_xlabel('number of Voronoi neighbours')
                        ax.set_ylabel('avg Voronoi cell surface area ($\AA^2$)')
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
                                ax_aggregate_inner_area.set_title('{condition} aggregate inner Leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp))
                                inner_species_counter += 1
                            else:
                                ax_aggregate_inner_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name])
                        else:
                            list_aggregate_frequency_data_outer_leaflet.extend(frequency_of_neighbour_count * [num_neighbours])
                            if outer_species_counter == 0:
                                ax_aggregate_outer_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name],label=molecular_species_name)
                                ax_aggregate_outer_area.set_title('{condition} aggregate outer Leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp))
                                outer_species_counter += 1
                            else:
                                ax_aggregate_outer_area.scatter(num_neighbours,average_surface_area,color=color_dict[molecular_species_name])
                                
                    ax2.bar(numpy.array(num_neighbours_list) - 0.3, frequency_of_neighbour_count_list)
                    if general_plots_xmax:
                        ax2.set_xlim(-1,general_plots_xmax)
                    else:
                        ax2.set_xlim(-1,12)
                    ax2.set_ylabel('frequency')
                    ax2.set_xlabel('number of Voronoi neighbours')
                    for axis in [ax,ax2]:
                        axis.set_title('{condition} {leaflet_flag} Leaflet {species} ({time_stamp})'.format(leaflet_flag = leaflet_name, species = molecular_species_name, condition = sim_name, time_stamp = time_stamp))
                    subplot_index += 1
                    #print 'subplot_index:', subplot_index
                    frequency_count_molecular_species += sum(frequency_of_neighbour_count_list)
                assert frequency_count_molecular_species == species_count_dictionary[molecular_species_name], "The neighbour frequency count for {mol_species} does not match the total molecules of this type in the system. Got {actual_count} instead.".format(mol_species = molecular_species_name, actual_count = frequency_count_molecular_species)
                
            for axis in [ax_aggregate_inner_area, ax_aggregate_outer_area]:
                axis.legend(loc=2,fontsize=10,scatterpoints=1)
                axis.set_xlabel('number of Voronoi neighbours')
                axis.set_ylabel('avg Voronoi cell surface area ($\AA^2$)')
                if general_plots_xmax:
                    axis.set_xlim(0,general_plots_xmax)
                else:
                    axis.set_xlim(0,14)
                axis.set_ylim(0,350)
                
            for axis in [ax_aggregate_inner_freq,ax_aggregate_outer_freq]:
                axis.set_xlabel('number of Voronoi neighbours')
                if general_plots_xmax:
                    axis.set_xlim(0,general_plots_xmax)
                else:
                    axis.set_xlim(0,14)
                axis.set_ylabel('Frequency')
                if aggregate_ymax:
                    axis.set_ylim(0,aggregate_ymax)
                else:
                    axis.set_ylim(0,2000)
                
            assert len(list_aggregate_frequency_data_inner_leaflet) + len(list_aggregate_frequency_data_outer_leaflet) == sum(species_count_dictionary.itervalues()), "The total number of Voronoi cells for which there are requency values should match the total number of particles in the system."

            bins = numpy.arange(14) + 0.1
            ax_aggregate_inner_freq.set_title('{condition} aggregate inner leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp))
            ax_aggregate_inner_freq.hist(list_aggregate_frequency_data_inner_leaflet,bins = bins)
            ax_aggregate_outer_freq.set_title('{condition} aggregate outer leaflet ({time_stamp})'.format(condition = sim_name, time_stamp = time_stamp))
            ax_aggregate_outer_freq.hist(list_aggregate_frequency_data_outer_leaflet,bins = bins)


            list_subplot_numbers_current_time_stamp = numpy.array(list_subplot_numbers_current_time_stamp) + 2
            aggregate_subplot_counter += 2








#            ax = self.fig.add_subplot(100,2,self.subplot_number)
#            ax.set_xlabel('num neighbours')
#            ax.set_ylabel('avg Voronoi cell\n surface area ($\AA^2$)')
#            ax.set_title(leaflet_name + ' leaflet ' + lipid_name)
#            self.subplot_number += 1
        self.fig.set_size_inches(fig_width_inches,fig_height_inches) 
        #self.fig.subplots_adjust(hspace = 0.7, wspace = 0.3)
        aggregate_figure_object.set_size_inches(fig_width_inches,10) 
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

def plot_sample_Voronoi_diagrams_zoom(matplotlib_figure_object,list_Voronoi_indices,dict_key_Voronoi_data,plot_title,dict_data,dengue_condition=None,debug_generators=None,generator_array=None):
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
    matplotlib_figure_object.savefig('sim35_outer_leaflet_pydata_london_2015.png', dpi = 300)

def area_per_molecule_plotting(figure_object,list_frame_numbers,list_percent_surface_area_reconstitution=None,list_percent_surface_area_reconstitution_inner_leaflet=None,protein_present=None,simulation_title=None,dictionary_headgroup_data=None,list_percent_surface_area_reconstitution_from_lipids_only=None,list_percent_surface_area_reconstitution_from_lipids_only_inner_leaflet=None,list_percent_surface_area_reconstitution_from_proteins_only=None,list_percent_surface_area_reconstitution_from_proteins_only_inner_leaflet=None,control_condition=None,dengue_condition=None):
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
