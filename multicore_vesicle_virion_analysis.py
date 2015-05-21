'''September 2012: General module for multi-core analysis of vesicle/virion trajectories using MDA and the multiprocessing module. Also have a specific interest in analysing lipid leaflet properties in the sim35/36 vesicles/virions. Among other things, I will probably want to be able to track the number of lipids of each type in a given leaflet and the area/lipid as the simulations progress. It might also be useful to be able to generate a (VMD ?) movie of some sort for visual verification of the assignment of species to leaflets.'''

import MDAnalysis, numpy, matplotlib, glob, time, multiprocessing
import cPickle as pickle
import cProfile
matplotlib.use('Agg') #deal with remote X11 issues

#because my large virion trajectories are continuously growing with various checkpoint parts and I also want to include the previous vesicle equilibrations in my analysis, I'd like to define a function that will return an MDA universe object that encompasses all of the trajectories of interest:
def generate_merged_MDA_universe(coordinate_file_path,ordered_list_of_trajectory_paths):
    '''Function to generate and return an MDAnalysis Universe object that encompasses a list of trajectories of interest for unified parsing. I suspect the list of trajectory paths must be properly ordered for the parsing to occur in a linear order.'''
    universe_object = MDAnalysis.Universe(coordinate_file_path, ordered_list_of_trajectory_paths)
    return universe_object

def produce_list_trajectories(path_string,glob_search_string):
    '''Function to produce a list of trajectories given a path string and a glob search string for the file names at that path. Include the final slash in the path_string'''
    list_of_trajectories = glob.glob(path_string + glob_search_string)
    return list_of_trajectories

def trajectory_list_order_checking(list_of_trajectories):
    '''A simple utility function for checking that trajectory lists used in the code are in the proper order (i.e., part_1.xtc, part_2.xtc, part_3.xtc) for loading into Universe objects.'''
    for trajectory in list_of_trajectories:
        print trajectory

def frame_iterator(universe_object):
    '''A simple test function to iterate through the frames of an MDAnalysis Universe object trajectory and print to the terminal along with the multiprocessing process name so that I can see that things are working properly.'''
    for ts in universe_object.trajectory[::]:
        print multiprocessing.current_process().name, 'frame : ', ts.frame

def track_PPCH_PO4_OD(queue_IPC,universe_object):
    '''Track the outer diameter of the vesicle/virion based on average PPCH PO4 particle radius from the COG of all lipids in the vesicle/virion (and just double that value in each frame to get the diameter instead of radius). Obviously this diameter measurement does not include the contribution of the large spike glycoproteins, but the idea is that I can track this particular parameter from the vesicle equilibration (where there are no proteins) right through the production runs of the virions to see how it changes before/after. This code also assumes that all PPCH residues remain in the outer leaflet, which remains to be tested rigorously (though PPCH looked a fair bit more stable in outer leaflet than some of the flip-flopping lipids early-on).'''
    import MDAnalysis.analysis.distances
    list_outer_diameter_values = [] #for storing the values calculated for this trajectory on this core
    PPCH_PO4_outer_leaflet_selection = universe_object.selectAtoms('resname PPCH and name PO4') #can set this before the traj loop as it will get dynamically updated during iteration
    #likewise for all lipid selection:
    all_lipid_selection = universe_object.selectAtoms('not resname RPO')
    for ts in universe_object.trajectory: 
        if ts.frame%100==0 or ts.frame ==1: #parse every 100th frame for now
            #select sphingomyelin phosphate particles as representative of the outer leaflet limits
            PPCH_PO4_outer_leaflet_coordinates = PPCH_PO4_outer_leaflet_selection.coordinates()
            PPCH_PO4_outer_leaflet_coordinates = numpy.array(PPCH_PO4_outer_leaflet_coordinates,dtype = 'float32')
            #select the COG of all lipids in the system as an overall centerpoint for measuring radii
            all_lipids_COG = all_lipid_selection.centerOfGeometry() #selection works because there's no solvent in sys
            all_lipids_COG = numpy.array(all_lipids_COG,dtype = 'float32')
            #determine all possible radial distances between PPCH phosphates and the centroid of the equilibrating vesicle / the virion
            PPCH_phosphate_radial_distances = MDAnalysis.analysis.distances.distance_array(all_lipids_COG[numpy.newaxis,:],PPCH_PO4_outer_leaflet_coordinates)
            #determine the average radial distance and double it; this should produce an estimate of the outer diameter in each frame:
            outer_diameter_estimate = (numpy.average(PPCH_phosphate_radial_distances) * 2) / 10.0 #also converting to nm units
            list_outer_diameter_values.append(outer_diameter_estimate) #so we'll have ordered list of OD values in nm
            print multiprocessing.current_process().name, 'frame : ', ts.frame #for tracking multi-core progress at terminal 
    #now, after the trajectory parsing has completed for this Universe object on this core, use the queue IPC object to send the numpy array of ordered OD values back to the parent process:
    array_OD_values = numpy.array(list_outer_diameter_values)
    #print 'array_OD_values.shape:',array_OD_values.shape #debugging
    queue_IPC.put(array_OD_values)

def capture_plot_PPCH_PO4_OD(Q1,Q2,Q3,Q4,Q5,Q6,from_pickle = 'no'):
    '''Capture the data from the child processes for the matching analysis function. Pickle the numpy array data to be safe and then plot the data with matplotlib.'''
    print 'Starting capture/plot code' #debugging
    #if not plotting from pickle files then generate numpy arrays directly from the multiprocessing queues:
    if from_pickle == 'no':
        #capture the numpy arrays sent back to queues from each core
        sim_33_array_PPCH_PO4_OD_values = Q1.get()
        sim_35_array_PPCH_PO4_OD_values = Q2.get()
        sim_36_array_PPCH_PO4_OD_values = Q3.get()
        sim_37_array_PPCH_PO4_OD_values = Q4.get()
        sim_38_array_PPCH_PO4_OD_values = Q5.get()
        sim_39_array_PPCH_PO4_OD_values = Q6.get()
        print 'starting to pickle' #debugging
        #pickle the numpy arrays to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim_33_array_PPCH_PO4_OD_values,open('sim33_PPCH_PO4_OD_tracking.p','wb'))
        pickle.dump(sim_35_array_PPCH_PO4_OD_values,open('sim35_PPCH_PO4_OD_tracking.p','wb'))
        pickle.dump(sim_36_array_PPCH_PO4_OD_values,open('sim36_PPCH_PO4_OD_tracking.p','wb'))
        pickle.dump(sim_37_array_PPCH_PO4_OD_values,open('sim37_PPCH_PO4_OD_tracking.p','wb'))
        pickle.dump(sim_38_array_PPCH_PO4_OD_values,open('sim38_PPCH_PO4_OD_tracking.p','wb'))
        pickle.dump(sim_39_array_PPCH_PO4_OD_values,open('sim39_PPCH_PO4_OD_tracking.p','wb'))
    elif from_pickle == 'yes': #for this case I simply retrieve the numpy arrays from the pickled files
        list_arrays= [pickle.load(open(x,'rb')) for x in ['sim33_PPCH_PO4_OD_tracking.p','sim35_PPCH_PO4_OD_tracking.p','sim36_PPCH_PO4_OD_tracking.p','sim37_PPCH_PO4_OD_tracking.p','sim38_PPCH_PO4_OD_tracking.p','sim39_PPCH_PO4_OD_tracking.p']]
        sim_33_array_PPCH_PO4_OD_values = list_arrays[0]
        sim_35_array_PPCH_PO4_OD_values = list_arrays[1]
        sim_36_array_PPCH_PO4_OD_values = list_arrays[2]
        sim_37_array_PPCH_PO4_OD_values = list_arrays[3]
        sim_38_array_PPCH_PO4_OD_values = list_arrays[4]
        sim_39_array_PPCH_PO4_OD_values = list_arrays[5]
    #now, regardless of which way the data is being read in by this function, I have a set of numpy arrays ready for matplotlib plotting:
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    #sim33 was the vesicle equilibration precursor to both sims 35/36, so produce merged arrays in each case:
    #sims33_35_combined_array = numpy.concatenate((sim_33_array_PPCH_PO4_OD_values[...,:31],sim_35_array_PPCH_PO4_OD_values))
    #sims33_36_combined_array = numpy.concatenate((sim_33_array_PPCH_PO4_OD_values[...,:31],sim_36_array_PPCH_PO4_OD_values))
    #for array in [sim_33_array_PPCH_PO4_OD_values,sim_35_array_PPCH_PO4_OD_values,sim_36_array_PPCH_PO4_OD_values]:
    #    matplotlib.pyplot.plot(numpy.arange(0,array.size)*10, array,linewidth = 2.0)
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('ytick', labelsize=8)
    matplotlib.pyplot.plot(numpy.arange(0,sim_33_array_PPCH_PO4_OD_values.size)*10/1000.0, sim_33_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#000000')
    matplotlib.pyplot.plot(numpy.arange(0,sim_35_array_PPCH_PO4_OD_values.size)*10/1000.0 + 0.3, sim_35_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#006600') #dark green
    matplotlib.pyplot.plot(numpy.arange(0,sim_36_array_PPCH_PO4_OD_values.size)*10/1000.0 + 0.3, sim_36_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#FF0000') #bright red
    matplotlib.pyplot.plot(numpy.arange(0,sim_37_array_PPCH_PO4_OD_values.size)*10/1000.0 + 0.3, sim_37_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#00FF00') #bright green (high T)
    matplotlib.pyplot.plot(numpy.arange(0,sim_38_array_PPCH_PO4_OD_values.size)*10/1000.0 + 0.3, sim_38_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#000066') #dark blue
    matplotlib.pyplot.plot(numpy.arange(0,sim_39_array_PPCH_PO4_OD_values.size)*10/1000.0 + 0.3, sim_39_array_PPCH_PO4_OD_values,linewidth = 4.0,c='#00CCFF') #light blue (high T)
    #matplotlib.pyplot.legend(('vesicle','virion (unrestrained)','virion (restrained)','virion (323K)','virion (FORS 295K)','virion (FORS 323K)'),loc = 0)
    matplotlib.pyplot.xlabel('Time ($\mu$s)',fontsize=14)
    matplotlib.pyplot.ylabel('PPCH ' + r'PO$_4$' + 'OD (nm)',fontsize=14)
    matplotlib.pyplot.subplots_adjust(left=0.15, right=0.9, top=0.90, bottom=0.15)
    matplotlib.pyplot.xlim(-0.2,5.3)
    matplotlib.pyplot.ylim(57,74)
    #matplotlib.pyplot.xlabel('Time (ns)')
    #matplotlib.pyplot.title('Tracking phosphate outer diameter \nbefore/after protein embedding to vesicles')
    #end_equilibration_time = sim_33_array_PPCH_PO4_OD_values.size * 10 #want to annotate the end point of the vesicle equil/start of virion
    #end_equilibration_OD = sim_33_array_PPCH_PO4_OD_values[-1]
    #matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(end_equilibration_time, end_equilibration_OD),xytext=(end_equilibration_time - 210, end_equilibration_OD + 4),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.05,frac=0.09))
    #fig.set_size_inches(11,7) #paper
    #fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/track_OD_extended/track_OD_extended_December_2013.png', dpi = 300)
    matplotlib.pyplot.subplots_adjust(left=0.20, right=0.97, top=0.90, bottom=0.15)
    fig.set_size_inches(4,4) 
    fig.savefig('/sansom/sc2/bioc1009/Dropbox/Tyler/flu_MS_figure_work/Fig_2_parts/Figure_2C_300_DPI.png', dpi = 300)

def track_multiple_protein_COGs_in_virion(queue_IPC,universe_object):
    '''Track the COG of each of the 107 membrane proteins in sim35 or sim36 and provide the results back to the queue via a list of COG values. Have to be conscious of the fact that I'm not deconvoluting the overall motion/rotation of the vesicles from the motion of their constituent proteins. This function was optimized using line_profiler and I removed a lot of internal looping and MDA stuff in favour of numpy slicing/operations.'''
    COG_values = []
    #generate all the index values for protein selections before iterating through trajectory
    #first the 80 index value ranges for HA:
    HA_atoms = 3621
    num_HA = 80
    list_index_ranges = [] #will store individual protein index values as tuples of boundaries here
    for index_range in range(num_HA): #generate 80 index ranges
        lower_index = index_range*HA_atoms
        upper_index = lower_index + HA_atoms
        list_index_ranges.append((lower_index,upper_index))
    #NA index values will start at the last upper index value:
    last_upper_index = list_index_ranges[-1][1]
    num_NA = 12
    NA_atoms = 4094
    for index_range in range(num_NA):
        lower_index = (index_range*NA_atoms) + last_upper_index
        upper_index = lower_index + NA_atoms
        list_index_ranges.append((lower_index,upper_index))
    #And similarly for M2:
    last_upper_index = list_index_ranges[-1][1]
    num_M2 = 15
    M2_atoms = 372
    for index_range in range(num_M2):
        lower_index = (index_range*M2_atoms) + last_upper_index
        upper_index = lower_index + M2_atoms
        list_index_ranges.append((lower_index,upper_index))
    #print 'final index:', list_index_ranges[-1][1],'(should be 344,388)'
    #so, should now be able to use this list_index_ranges to slice the numpy array of protein coordinates in each frame and calculate the 107 COG values (1 for each protein)
    #select all protein atoms in the system based on the index:
    protein_atom_selection = universe_object.selectAtoms('bynum 1:344388') #slight performance improvement if this is outside traj loop
    for ts in universe_object.trajectory[::100]: #parse every N frames for now
        #produce a numpy array of all the protein coordinates in the system:
        protein_coord_array = protein_atom_selection.coordinates() #should have 344388 rows and 3 columns
        #now iterate through the list of protein index ranges and slice the numpy array of coords to calculate the respective COG values:
        for lower_index,upper_index in list_index_ranges:
            #calculate COG of CG proteins as the average of all CG particle coordinates
            protein_COG = numpy.average(protein_coord_array[lower_index:upper_index,:],axis=0)
            #print protein_COG
            COG_values.append(protein_COG)
        print multiprocessing.current_process().name, 'frame : ', ts.frame #for tracking multi-core progress at terminal
    queue_IPC.put(COG_values)

def capture_plot_multiple_protein_COGs_in_virion(Q1,Q2,from_pickle = 'no'):
    '''Plot the data tracking the 107 membrane protein COGs in the virions, whether it be 'live' from the parallel MDA calculation or retrieved from pickled data. Sensitive to matplotlib version (3D plotting stuff)--using nanotube rather than my machine at the moment for this reason.'''
    if from_pickle == 'no': #first run-through of calculation; grab the data directly from the queue objects
        sim_35_COG_list = Q1.get()
        sim_36_COG_list = Q2.get()
        #and pickle the data to be safe and save time later:
        pickle.dump(sim_35_COG_list,open('sim35_COG_list.p','wb'))
        pickle.dump(sim_36_COG_list,open('sim36_COG_list.p','wb'))
    elif from_pickle == 'yes':
        #in this case retrieve the list of COG arrays from the pickled data:
        sim_35_COG_list = pickle.load(open('sim35_COG_list.p','rb'))
        sim_36_COG_list = pickle.load(open('sim36_COG_list.p','rb'))
    #now, regardless of which method was used to retrieve data, there are two lists of COG numpy arrays ready for plotting
    sim_35_COG_array = numpy.array(sim_35_COG_list)
    sim_36_COG_array = numpy.array(sim_36_COG_list)
#    print sim_35_COG_array.shape
#    print 'Expecting around 13,000 at current progression.'
    print 'Arrays available for plotting' #debugging
    import matplotlib.pyplot
    from mpl_toolkits.mplot3d import Axes3D
    matplotlib.rcParams.update({'font.size': 22})
    fig = matplotlib.pyplot.figure()
    #plot the unrestrained data (sim#35):
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_xticks(numpy.arange(400,1200,200))
    ax1.set_yticks(numpy.arange(400,1200,200))
    ax1.set_zticks(numpy.arange(0,700,200))
    ax1.scatter(sim_35_COG_array[:,0].ravel(),sim_35_COG_array[:,1].ravel(),sim_35_COG_array[:,2].ravel(),zdir='z',c='r',marker='d')
#    ax1.set_title('Unrestrained non-Forss\n virion protein COG\n motion (1251 ns)')
    ax1.set_xlabel('\n\nx' + r'($\AA$)')
    ax1.set_ylabel('\n\ny' + r'($\AA$)')
    ax1.set_zlabel('\n\nz' + r'($\AA$)')
    #plot the restrained data:
    fig2 = matplotlib.pyplot.figure()
    ax2 = fig2.add_subplot(111, projection = '3d')
    ax2.set_xticks(numpy.arange(400,1200,200))
    ax2.set_yticks(numpy.arange(400,1200,200))
    ax2.set_zticks(numpy.arange(0,700,200))
    ax2.scatter(sim_36_COG_array[:,0].ravel(),sim_36_COG_array[:,1].ravel(),sim_36_COG_array[:,2].ravel(),zdir='z',c='r',marker='d')
#    ax2.set_title('Pull-code-restrained\n non-Forss virion protein COG\n motion (880 ns)')
    ax2.set_xlabel('\n\nx' + r'($\AA$)')
    ax2.set_ylabel('\n\ny' + r'($\AA$)')
    ax2.set_zlabel('\n\nz' + r'($\AA$)')
    fig.set_size_inches(10,10)
    fig2.set_size_inches(10,10)
    fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/COG_motion_sims35_36/sim35_COG_motion_extended_Jan_2014.png',dpi=300)
    fig2.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/COG_motion_sims35_36/sim36_COG_motion_extended_Jan_2014.png',dpi=300)
    
def leaflet_asymmetry_analysis(queue_IPC,universe_object):
    '''Code for analyzing various leaflet asymmetry properties in the vesicle/virion simulations.'''
    #start off by defining the selection string to be used for leaflet networking in MDA/networkx:
    #I'm going to test excluding CHOL from the networking process because I suspect it is bridging networks with its burial position:
    leaflet_selection_string = 'name PO4' #should cover the lipids except CHOL (but would need adjusting for Forssman)
    import MDAnalysis.analysis.leaflet
    list_num_leaflets = [] #for tracking the number of leaflets identified by LeafletFinder in each frame using my current cutoff value
    #iterate through the trajectory and try to establish leaflet assignments for lipids:
    for ts in universe_object.trajectory[::10]: #every N frames for testing
        L = MDAnalysis.analysis.leaflet.LeafletFinder(universe_object, leaflet_selection_string, cutoff=19.0)
        list_leaflet_groups = L.groups()
        num_identified_leaflets = len(list_leaflet_groups)
        #now print out the number of leaflets identified (of course, the target is to have 2--no more, no less--in each frame)
        #also append the number of identified leaflets to a list which can eventually be tracked in a plot:
        list_num_leaflets.append(num_identified_leaflets)
        print 'Number of identified leaflets:', num_identified_leaflets
        print multiprocessing.current_process().name, 'frame : ', ts.frame #for tracking multi-core progress at terminal
    #when this core has finished parsing this trajectory, send the list of num_identified_leaflets back to the queue object:
    queue_IPC.put(list_num_leaflets)

def capture_plot_leaflet_asymmetry_analysis(Q1,Q2,Q3,from_pickle = 'no'):
    if from_pickle == 'no': #retrieving the live data from the queue objects
        sim33_num_leaflets_list = Q1.get()
        sim35_num_leaflets_list = Q2.get()
        sim36_num_leaflets_list = Q3.get()
        sim33_array_num_leaflets = numpy.array(sim33_num_leaflets_list)
        sim35_array_num_leaflets = numpy.array(sim35_num_leaflets_list)
        sim36_array_num_leaflets = numpy.array(sim36_num_leaflets_list)
        #pickle the data for safety/future use:
        pickle.dump(sim33_array_num_leaflets,open('sim33_leaflet_count_tracking.p','wb'))
        pickle.dump(sim35_array_num_leaflets,open('sim35_leaflet_count_tracking.p','wb'))
        pickle.dump(sim36_array_num_leaflets,open('sim36_leaflet_count_tracking.p','wb'))
    elif from_pickle == 'yes': #retrieving the data from pickled files
        sim33_array_num_leaflets = pickle.load(open('sim33_leaflet_count_tracking.p','rb'))
        sim35_array_num_leaflets = pickle.load(open('sim35_leaflet_count_tracking.p','rb'))
        sim36_array_num_leaflets = pickle.load(open('sim36_leaflet_count_tracking.p','rb'))
    #now regardless of the method used to retrieve data, numpy arrays are available for matplotlib plotting
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    #concatenate the vesicle equilibration leaflet number tracking with the respective virion simulations for a complete plot start/finish
    sim35_merged_array = numpy.concatenate((sim33_array_num_leaflets,sim35_array_num_leaflets))
    sim36_merged_array = numpy.concatenate((sim33_array_num_leaflets,sim36_array_num_leaflets))
    for array in [sim35_merged_array,sim36_merged_array]: 
        matplotlib.pyplot.plot(numpy.arange(0,array.size),array, linewidth = 2.0)
    matplotlib.pyplot.ylim((0,6))
    matplotlib.pyplot.ylabel('Number of leaflets assigned\n' +r'(PO$_4$ 19$\AA$ cutoff)')
    matplotlib.pyplot.xlabel('Time (ns)')
    matplotlib.pyplot.title('Assessing confidence in LeafletFinder \nbased on number of assigned leaflets')
    matplotlib.pyplot.legend(('unrestrained (sims 33_35)','restrained (sims 33_36)'),loc=0)
    fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/leaflet_analysis/track_num_assigned_leaflets.png',dpi=300)

def leaflet_asymmetry_analysis_2(queue_IPC, universe_object):
    '''Version 2: Code for analyzing various leaflet asymmetry properties in the vesicle/virion simulations.'''
    #be careful with selections: remember that my intention is to parse the vesicle before AND after protein incorporation--this means the number of lipids changes and index numbers change drastically during the overall time course
    #produce useful atom selections outside of the trajectory loop where possible:
    #start off by defining the selection string to be used for leaflet networking in MDA/networkx:
    #I'm going to test excluding CHOL from the networking process because I suspect it is bridging networks with its burial position:
    leaflet_selection_string = 'name PO4' #should cover the lipids except CHOL (but would need adjusting for Forssman)
    #the total number of lipids in the system will vary as it progresses through equilibration/protein incorporation, so produce a sel:
    #total_lipid_excl_chol = universe_object.selectAtoms('not resname RPO and not resname CHOL') #this sel works because no solvent in trajs
    #testing an alternate selection:
    total_lipid_excl_chol = universe_object.selectAtoms('resname PPCH or resname POPS or resname DOPX or resname DOPE') #this selection works, the previous one didn't
    total_lipid_incl_chol = universe_object.selectAtoms('not resname RPO') #will probably use this to determine vesicle/virion centroid
    import MDAnalysis.analysis.leaflet, MDAnalysis.analysis.distances
    #initialize a set of lists for tracking parameters to be plotted after parsing trajectories:
    list_num_leaflets = [] #for tracking the number of leaflets identified by optimize_cutoff in each frame
    list_total_lipids_per_leaflet = [] #will append a 3-element tuple in each frame corresponding to the % of total lipids found in the 3 most populous leaflets (simply place 0 if one or two if the leaflets is not assigned)
    list_phosphate_dist_vesicle_centroid_per_leaflet = [] #will use similar 3-element tuples in each frame to track the average phosphate distance from the vesicle centroid in each of the 3 most populous leaflets
    list_optimized_cutoff_values = [] #simply append the optimized cutoff value from each frame to this list
    #iterate through the trajectory and do required work:
    for ts in universe_object.trajectory[::10]: #every N frames for testing
        #determine the optimal cutoff for network propagation using the built-in MDA function; will have to use a pretty wide range of values as the vesicle starts off very large and then shrinks down before protein incorporation and additional simulation:
        try:
            optimal_cutoff, num_leaflets = MDAnalysis.analysis.leaflet.optimize_cutoff(universe_object,leaflet_selection_string,dmin=10.0,dmax=40.0, step=1.0,max_imbalance = 0.5) #note that I'm inflating max_imbalance because of leaflet size discrepencies in vesicle
        except:
            optimal_cutoff, num_leaflets = (0,0) #in case of problem with algorithm, which seems to occur for sim33
        #append the optimized cutoff and the number of leaflets to their respective lists:
        list_num_leaflets.append(num_leaflets)
        list_optimized_cutoff_values.append(optimal_cutoff)
        #now use LeafletFinder() to parse out the leaflet assignments using the optimized cutoff
        L = MDAnalysis.analysis.leaflet.LeafletFinder(universe_object, leaflet_selection_string, cutoff=optimal_cutoff)
        #determine the total number of non-CHOL lipid residues in the system:
        total_non_CHOL_lipid_residues = total_lipid_excl_chol.numberOfResidues()
        print 'total non CHOL lipid residues:', total_non_CHOL_lipid_residues
        #produce a dictionary containing the number of residues in each leaflet:
        leaflet_size_dict = L.sizes()
        #initialize a list for appending the % total non-CHOL lipid in each leaflet:
        percent_total_lipid_non_chol_sublist = [] #this will, in turn, be appended as a nested list to the overall list
        #convert the number of residues in each leaflet to a percentage of the total available (non-CHOL) residues in the system:
        for num_residues in leaflet_size_dict.values():
            percentage_total_residues_per_leaflet = (float(num_residues) / float(total_non_CHOL_lipid_residues)) * 100.0
            percent_total_lipid_non_chol_sublist.append(percentage_total_residues_per_leaflet)
            print 'Percentage in current leaflet:', percentage_total_residues_per_leaflet #debugging
        #now filter the sublist to only contain 3 entries--1 each for the first (largest) 3 leaflets
        if len(percent_total_lipid_non_chol_sublist) == 3: #what we want
            pass #don't need to do anything here
        elif len(percent_total_lipid_non_chol_sublist) < 3: #pad the empty leaflet entries with 0
            while len(percent_total_lipid_non_chol_sublist) < 3:
                percent_total_lipid_non_chol_sublist.append(0) #0% non-chol lipids present in non-existent leaflet
        elif len(percent_total_lipid_non_chol_sublist) > 3: #slice it down to data for the 3 most populated leaflets
            percent_total_lipid_non_chol_sublist = percent_total_lipid_non_chol_sublist[0:3]
        #finally append the 3 % values to the overall list:
        list_total_lipids_per_leaflet.append(percent_total_lipid_non_chol_sublist) #will have N rows by 3 columns (1 column for each leaflet/frame)
        #determine the centroid of all lipids in the vesicle in this frame:
        overall_centroid_coordinates = total_lipid_incl_chol.centerOfGeometry()
        overall_centroid_coordinates = numpy.array(overall_centroid_coordinates, dtype = 'float32')
        #create a sublist for appending the avg PO4 radial distance, which will then be nested into an overall list:
        avg_PO4_radial_dist_sublist = []
        #produce a numpy array of coordinates for all PO4 particles in each leaflet in this frame:
        for atom_group in L.groups():
            try:
                array_leaflet_PO4_coordinates = atom_group.PO4.positions
            except: #presumably the exception results from having a single PO4 atom in a 'leaflet'
                array_leaflet_PO4_coordinates = atom_group.PO4.position
                array_leaflet_PO4_coordinates = array_leaflet_PO4_coordinates[numpy.newaxis,:]
            array_leaflet_PO4_coordinates = numpy.array(array_leaflet_PO4_coordinates,dtype='float32')
            #now generate an array of distances between each PO4 particle in current leaflet and overall vesicle centroid:
            dist_array = MDAnalysis.analysis.distances.distance_array(overall_centroid_coordinates[numpy.newaxis,:],array_leaflet_PO4_coordinates)
            #now calculate the average PO4 distance (r) from the current vesicle COG for this leaflet in this frame:
            avg_PO4_radial_distance_this_leaflet = numpy.average(dist_array)/10.0
            avg_PO4_radial_dist_sublist.append(avg_PO4_radial_distance_this_leaflet)
            print 'Current leaflet average PO4 radius from vesicle COG (nm):', avg_PO4_radial_distance_this_leaflet
        #now filter the sublist to only contain 3 entries--1 each for the first (largest) 3 leaflets
        if len(avg_PO4_radial_dist_sublist) == 3: #what we want
            pass
        elif len(avg_PO4_radial_dist_sublist) < 3: #pad empty leaflet entries with 0
            while len(avg_PO4_radial_dist_sublist) < 3:
                avg_PO4_radial_dist_sublist.append(0)
        elif len(avg_PO4_radial_dist_sublist) > 3: #slice it down to top 3 entries in this case
            avg_PO4_radial_dist_sublist = avg_PO4_radial_dist_sublist[0:3]
        #finally append the PO4 radial distance sublist to the overall list:
        list_phosphate_dist_vesicle_centroid_per_leaflet.append(avg_PO4_radial_dist_sublist)

        print 'Number of leaflets:', num_leaflets #debugging
        print 'Dict of leaflet sizes:', leaflet_size_dict #debugging
        print 'Optimized cutoff value (A):', optimal_cutoff #debugging
        print multiprocessing.current_process().name, 'frame : ', ts.frame #for tracking multi-core progress at terminal
    #after the trajectory iteration completes, put the 4 lists to be used for the 4 plots into the queue object:
    queue_IPC.put([list_num_leaflets,list_optimized_cutoff_values,list_total_lipids_per_leaflet,list_phosphate_dist_vesicle_centroid_per_leaflet])

def capture_plot_leaflet_asymmetry_analysis_2(Q1,Q2,Q3,from_pickle = 'no'):
    if from_pickle == 'no': #retrieve data from queue objects
        #data lists should have format: [list_num_leaflets,list_optimized_cutoff_values,list_total_lipids_per_leaflet,list_phosphate_dist_vesicle_centroid_per_leaflet]
        sim_33_data_list = Q1.get() 
        sim_35_data_list = Q2.get() 
        sim_36_data_list = Q3.get() 
        #pickle the data for safety/future use:
        pickle.dump(sim_33_data_list,open('sim_33_leaflet_analysis_2.p','wb'))
        pickle.dump(sim_35_data_list,open('sim_35_leaflet_analysis_2.p','wb'))
        pickle.dump(sim_36_data_list,open('sim_36_leaflet_analysis_2.p','wb'))
    elif from_pickle == 'yes': #retrieve data from pickled objects
        sim_33_data_list = pickle.load(open('sim_33_leaflet_analysis_2.p','rb'))
        sim_35_data_list = pickle.load(open('sim_35_leaflet_analysis_2.p','rb'))
        sim_36_data_list = pickle.load(open('sim_36_leaflet_analysis_2.p','rb'))
    #now, regardless of which method was used to retrieve the lists of data, they are avaible for manipulation/plotting
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    #unpack the nested lists for sim33:
    sim33_list_num_leaflets,sim33_list_optimized_cutoff_values,sim33_list_total_lipids_per_leaflet,sim33_list_phosphate_dist_vesicle_centroid_per_leaflet = (element for element in sim_33_data_list)
    #likewise for sims35/36:
    sim35_list_num_leaflets,sim35_list_optimized_cutoff_values,sim35_list_total_lipids_per_leaflet,sim35_list_phosphate_dist_vesicle_centroid_per_leaflet = (element for element in sim_35_data_list)
    sim36_list_num_leaflets,sim36_list_optimized_cutoff_values,sim36_list_total_lipids_per_leaflet,sim36_list_phosphate_dist_vesicle_centroid_per_leaflet = (element for element in sim_36_data_list)
    #produce merged lists for sims33/35 and sims33/36:
    sims33_35_list_num_leaflets = sim33_list_num_leaflets + sim35_list_num_leaflets
    sims33_35_list_optimized_cutoff_values = sim33_list_optimized_cutoff_values + sim35_list_optimized_cutoff_values
    sims33_35_list_total_lipids_per_leaflet = sim33_list_total_lipids_per_leaflet + sim35_list_total_lipids_per_leaflet
    sims33_35_list_phosphate_dist_vesicle_centroid_per_leaflet = sim33_list_phosphate_dist_vesicle_centroid_per_leaflet + sim35_list_phosphate_dist_vesicle_centroid_per_leaflet
    sims33_36_list_num_leaflets = sim33_list_num_leaflets + sim36_list_num_leaflets
    sims33_36_list_optimized_cutoff_values = sim33_list_optimized_cutoff_values + sim36_list_optimized_cutoff_values
    sims33_36_list_total_lipids_per_leaflet = sim33_list_total_lipids_per_leaflet + sim36_list_total_lipids_per_leaflet
    sims33_36_list_phosphate_dist_vesicle_centroid_per_leaflet = sim33_list_phosphate_dist_vesicle_centroid_per_leaflet + sim36_list_phosphate_dist_vesicle_centroid_per_leaflet

    matplotlib.pyplot.subplots_adjust(left=0.15, right=0.9, top=0.95, bottom=0.1)

    #start setting up the plots:
    ax1 = fig.add_subplot(421) #track num leaflets for sims33/35
    ax1.plot(numpy.arange(0,numpy.array(sims33_35_list_num_leaflets).size),numpy.array(sims33_35_list_num_leaflets),linewidth = 1.0)
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('# Assigned leaflets\nexcluding CHOL')
    ax1.set_ylim((0,6))
    ax1.set_xlim((0,1100))
    ax1.set_title('Sims 33/35')
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 2),xytext=(300, 0.4),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))
    
    ax2 = fig.add_subplot(422) #track num leaflets for sims33/36
    ax2.plot(numpy.arange(0,numpy.array(sims33_36_list_num_leaflets).size),numpy.array(sims33_36_list_num_leaflets),linewidth = 1.0)
    ax2.set_xlabel('Time (ns)')
    #ax2.set_ylabel('# Assigned leaflets\nexcluding CHOL')
    ax2.set_ylim((0,6))
    ax2.set_xlim((0,1100))
    ax2.set_title('Sims 33/36')
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 2),xytext=(300, 0.4),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))

    ax3 = fig.add_subplot(423) #track optimized cutoff for sims33/35
    ax3.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),sims33_35_list_optimized_cutoff_values,linewidth = 1.0)
    ax3.set_xlabel('Time (ns)')    
    ax3.set_ylabel('optimized\nnetworking\ncutoff ' + r'($\AA$)')
    ax3.set_xlim((0,1100))
    ax3.set_ylim((0,40))
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 15),xytext=(300, 4),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))

    ax4 = fig.add_subplot(424) #track optimized cutoff for sim33/36
    ax4.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),sims33_36_list_optimized_cutoff_values,linewidth = 1.0)
    ax4.set_xlabel('Time (ns)')    
    ax4.set_xlim((0,1100))
    ax4.set_ylim((0,40))
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 15),xytext=(300, 4),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))
    
    ax5 = fig.add_subplot(425) #track %total non-CHOL lipids in the 3 most populous assigned leaflets for sim33/35
    #plot data for leaflets 1,2,3 in order
    ax5.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_total_lipids_per_leaflet)[:,0].ravel(),linewidth = 1.0, label = 'leaflet 1')
    ax5.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_total_lipids_per_leaflet)[:,1].ravel(),linewidth = 1.0, label = 'leaflet 2')
    ax5.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_total_lipids_per_leaflet)[:,2].ravel(),linewidth = 1.0, label = 'leaflet 3')
    ax5.set_xlabel('Time (ns)')    
    ax5.set_ylabel('% Total non-CHOL lipid')
    ax5.legend(loc=0)
    ax5.set_xlim((0,1100))
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 50),xytext=(300, 30),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))

    ax6 = fig.add_subplot(426) #same as 425 but for sim33/36
    #plot data for leaflets 1,2,3 in order
    ax6.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_total_lipids_per_leaflet)[:,0].ravel(),linewidth = 1.0, label = 'leaflet 1')
    ax6.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_total_lipids_per_leaflet)[:,1].ravel(),linewidth = 1.0, label = 'leaflet 2')
    ax6.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_total_lipids_per_leaflet)[:,2].ravel(),linewidth = 1.0, label = 'leaflet 3')
    ax6.set_xlabel('Time (ns)')    
    ax6.legend(loc=0)
    ax6.set_xlim((0,1100))
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 50),xytext=(300, 30),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))

    ax7 = fig.add_subplot(427) #track the average PO4 radial distance from virion centroid for each of the top 3 most populous leaflets (sims33/35):
    #once again, plot data for leaflets 1,2,3 in order
    ax7.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,0].ravel(),linewidth = 1.0, label = 'leaflet 1',zorder=3)
    ax7.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,1].ravel(),linewidth = 1.0, label = 'leaflet 2',zorder=2)
    ax7.plot(numpy.arange(0,numpy.array(sims33_35_list_optimized_cutoff_values).size),numpy.array(sims33_35_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,2].ravel(),linewidth = 1.0, label = 'leaflet 3',zorder=1)
    ax7.set_xlabel('Time (ns)')
    ax7.legend(loc=0)
    ax7.set_xlim((0,1100))
    ax7.set_ylabel('Average PO4 radius\nfrom vesicle centroid (nm)')
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 25),xytext=(400, 16),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))

    ax8 = fig.add_subplot(428) #similar to 427 but for sim33/36
    #once again, plot data for leaflets 1,2,3 in order
    ax8.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,0].ravel(),linewidth = 1.0, label = 'leaflet 1',zorder=3)
    ax8.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,1].ravel(),linewidth = 1.0, label = 'leaflet 2',zorder =2)
    ax8.plot(numpy.arange(0,numpy.array(sims33_36_list_optimized_cutoff_values).size),numpy.array(sims33_36_list_phosphate_dist_vesicle_centroid_per_leaflet)[:,2].ravel(),linewidth = 1.0, label = 'leaflet 3',zorder=1)
    ax8.set_xlabel('Time (ns)')
    ax8.legend(loc=0)
    ax8.set_xlim((0,1100))
    matplotlib.pyplot.annotate('end vesicle equil;\nadd proteins',xy=(300, 25),xytext=(400, 16),arrowprops=dict(facecolor='black', width=2.0,shrink = 0.10,frac=0.15))



    fig.set_size_inches(8.5,16)
    fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/leaflet_analysis/leaflet_analysis_2.png',dpi = 300)


def stratification_tracking(queue_IPC,universe_object):
    '''Parallel version of the serial-optimized code. This function is intended as an alternative to an MDA LeafletFinder()-based approach that allows for tracking of the probability (via contour plotting) of finding a member of given lipid species at a given distance (layer/strata) away from the centroid of the entire vesicle (excluding RPO/proteins). I'm hoping this may also perform much more rapidly than the painfully slow code required for network optimization and LeafletFinder(). It will be especially important to verify that PPCH remains almost entirely in the outer leaflet as I have been using PPCH phosphates to track the outer diameter of the lipid component of the virions.'''
    import MDAnalysis.analysis.distances
    #regardless of whether I'm parsing sim33 or sims35/36, I can use the same lipid residue selection for assessing the overall COG of the system (excluding RPO particles and protein particles):
    #in preparation for Biophysics 2013, I'm going to try including FORS here as well for sim38:
    overall_lipid_selection_string = 'resname POPS or resname DOPX or resname DOPE or resname CHOL or resname PPCH or resname FORS'
    overall_lipid_selection = universe_object.selectAtoms(overall_lipid_selection_string)
    #now, as I'm going to focus on either phosphate or ROH (CHOL) particles for these distance measurements, I can define these before the loop (and because I want a separate contour plot for each lipid species):
    #add / select only those currently needed for Biophysics 2013:
    POPS_PO4_selection = universe_object.selectAtoms('resname POPS and name PO4')
    #DOPX_PO4_selection = universe_object.selectAtoms('resname DOPX and name PO4')
    #DOPE_PO4_selection = universe_object.selectAtoms('resname DOPE and name PO4')
    CHOL_ROH_selection = universe_object.selectAtoms('resname CHOL and name ROH')
    PPCH_PO4_selection = universe_object.selectAtoms('resname PPCH and name PO4')
    #will almost certainly want to include FORS for tracking as well, at some point, but for now just use the lipid species that I can compare across many conditions
    #now, define a function that, in each frame, will return matplotlib-ready binned positional probabilities for the set of phosphates/ROH for a given lipid species:
    def calc_positional_probabilities(MDA_selection,bin_step,overall_COG):
        '''the bin_step (in nm) can be reduced (smaller histogram bins for particle distances from COG) to increase resolution'''
        #I'm going to want the set of PO4 or ROH coords for a given MDAnalysis lipid selection of interest:
        coordinate_array = MDA_selection.coordinates() #an Nx3 numpy array of coordinates
        #now produce an array of the distances from the above coordinates to the overall lipid COG, and convert to nm units
        array_of_distances = MDAnalysis.analysis.distances.distance_array(numpy.array(overall_COG[numpy.newaxis,:],dtype='float32'),numpy.array(coordinate_array,dtype='float32'))/10.0
        #I think the above array will have shape (1,N) or size N; regardless the array of distances will be flattened when 'histogrammed'
        #now use the built-in numpy histogram method to deal with the binning stuff:
        array_bin_counts, array_bin_edges = numpy.histogram(array_of_distances,bins=numpy.arange(0,80,bin_step))
        #now convert the array of bin counts to percentages:
        array_bin_percentages = ( array_bin_counts / float(numpy.sum(array_bin_counts)) ) * 100.0
        x_array = numpy.array([ts.frame]) #just 1 'band' of the final contour plot per frame
        z_array = array_bin_percentages 
        midpoint_values = []
        for index in range(array_bin_edges.size-1):
            current_midpoint_value = (array_bin_edges.flat[index]+ array_bin_edges.flat[index+1]) /2.0
            midpoint_values.append(current_midpoint_value)
        y_array = numpy.array(midpoint_values)
        #now return these arrays in a tuple
        return (x_array,y_array,z_array)

    #and now I will want to initialize some lists for storing the numpy arrays of data from each parsed frame for a given lipid species:
    POPS_data = [] #this and subsequent species data arrays will have format x_array,y_array,z_array repeated for each time frame (x value)
    #DOPX_data = []
    #DOPE_data = []
    CHOL_data = []
    PPCH_data = []
        
    for ts in universe_object.trajectory: 
        #I will definitely need the overall COG of the lipids in the sys in each parsed frame:
        overall_lipid_COG = overall_lipid_selection.centerOfGeometry() #a 1x3 numpy array for the overall lipid COG coords
        #determine the histogrammed distance probability data formatted for matplotlib for each of the lipids of interest in this frame:
        index = 0
        #for MDA_selection in [POPS_PO4_selection,DOPX_PO4_selection,DOPE_PO4_selection,CHOL_ROH_selection,PPCH_PO4_selection]:
        for MDA_selection in [POPS_PO4_selection,CHOL_ROH_selection,PPCH_PO4_selection]:
            tuple_xyz_arrays = calc_positional_probabilities(MDA_selection,bin_step=0.5,overall_COG=overall_lipid_COG)
            #append the arrays to the appropriate data list for storage and eventual matplotlib plotting
            if index == 0:
                data_list = POPS_data
            #elif index == 1:
            #    data_list = DOPX_data
            #elif index == 2: 
            #    data_list = DOPE_data
            elif index == 1:
                data_list = CHOL_data
            elif index == 2: 
                data_list = PPCH_data
            index +=1
            for array in tuple_xyz_arrays:
                data_list.append(array) #so all of the data lists should always have the format x,y,z,x,y,z,etc...
        print multiprocessing.current_process().name, 'frame : ', ts.frame #for tracking multi-core progress at terminal
    #put the list of data back into the queue from this child process
    #queue_IPC.put([POPS_data,DOPX_data,DOPE_data,CHOL_data,PPCH_data])
    queue_IPC.put([POPS_data,CHOL_data,PPCH_data])

def capture_plot_stratification_tracking(Q1,Q2,Q3,Q4,Q5,Q6,from_pickle = 'no'):
    '''Updated to capture parallel data. Capture the stratification tracking data and produce appropriate contour plots using matplotlib.'''
    #Making adjustments for Biophysics 2013 panels:
    if from_pickle == 'no': # in this case, pickle the data for storage/future use
        sim33_list_of_data = Q1.get()
        sim35_list_of_data = Q2.get()
        sim36_list_of_data = Q3.get()
        sim37_list_of_data = Q4.get()
        sim38_list_of_data = Q5.get()
        sim39_list_of_data = Q6.get()
        pickle.dump(sim33_list_of_data,open('sim33_parallel_stratification_tracking.p','wb'))
        pickle.dump(sim35_list_of_data,open('sim35_parallel_stratification_tracking.p','wb'))
        pickle.dump(sim36_list_of_data,open('sim36_parallel_stratification_tracking.p','wb'))
        pickle.dump(sim37_list_of_data,open('sim37_parallel_stratification_tracking.p','wb'))
        pickle.dump(sim38_list_of_data,open('sim38_parallel_stratification_tracking.p','wb'))
        pickle.dump(sim39_list_of_data,open('sim39_parallel_stratification_tracking.p','wb'))
    elif from_pickle == 'yes': # in this case, retrieve the data from the pickled files
        sim33_list_of_data = pickle.load(open('sim33_parallel_stratification_tracking.p','rb'))
        sim35_list_of_data = pickle.load(open('sim35_parallel_stratification_tracking.p','rb'))
        sim36_list_of_data = pickle.load(open('sim36_parallel_stratification_tracking.p','rb'))
        sim37_list_of_data = pickle.load(open('sim37_parallel_stratification_tracking.p','rb'))
        sim38_list_of_data = pickle.load(open('sim38_parallel_stratification_tracking.p','rb'))
        sim39_list_of_data = pickle.load(open('sim39_parallel_stratification_tracking.p','rb'))
    '''
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.05)
    #iterate through the per-frame data and plot it for each lipid type in each simulation condition
    #also want to merge sim33 data into sim35/36 data to have continuous plots tracking from equilibration through to production
    zipped_sim33_35_list_of_data = zip(sim33_list_of_data,sim35_list_of_data)
    zipped_sim33_36_list_of_data = zip(sim33_list_of_data,sim36_list_of_data)
    #the above zipped lists should have formats liks this [(POPS_sim33_data, POPS_sim35_data), (DOPX_sim33_data, DOPX_sim35_data), ...]
    fig.set_size_inches(14,14)
    #plot sim33/35 data in the left column:
    ax_positions = [521, 523, 525, 527, 529]
    title_list = ['sim33/35\n POPS PO4 position','DOPX PO4 position','DOPE PO4 position','CHOL ROH position','PPCH PO4 position']
    for sim33_lipid_data_list, sim35_lipid_data_list in zipped_sim33_35_list_of_data:
        sim33_all_frame_numbers_array = numpy.array([array[0] for array in sim33_lipid_data_list[0::3]]) #every 3rd data value starting from first index is a frame number        
        sim35_all_frame_numbers_array = numpy.array([array[0] for array in sim35_lipid_data_list[0::3]]) 
        y_values_array = sim33_lipid_data_list[1] #these are all the same so can just grab one
        sim33_z_values_array = numpy.transpose(numpy.array(sim33_lipid_data_list[2::3])) #a bit of manipulation to get the shape of this array ready for matplotlib
        sim35_z_values_array = numpy.transpose(numpy.array(sim35_lipid_data_list[2::3])) 
        #combine the Z arrays so that sim33 values precede sim35 values
        z_values_array = numpy.concatenate((sim33_z_values_array,sim35_z_values_array),axis=1)
        #combine the frame number arrays after adding in the values for sim35 to avoid overlap with sim33:
        sim35_all_frame_numbers_array += sim33_all_frame_numbers_array[-1]
        x_values_array = numpy.concatenate((sim33_all_frame_numbers_array,sim35_all_frame_numbers_array))/ 10.0
        ax = fig.add_subplot(ax_positions.pop(0)) #position plots for each lipid type as you iterate through
        contour_plot = ax.contourf(x_values_array,y_values_array,z_values_array, levels = numpy.logspace(-8.0,6.0,num=50,base=2.0))
        cbar = matplotlib.pyplot.colorbar(contour_plot)
        cbar.ax.set_ylabel('% particles')
        ax.set_title(title_list.pop(0))
        ax.set_ylabel('Distance from\nlipid centroid (nm)')
        ax.set_ylim((10,45))
        ax.set_xlim((0,1100))
    ax.set_xlabel('Time (ns)')
    #and now likewise for sim33/36 data in the right column:
    ax_positions = [(5,2,2),(5,2,4),(5,2,6),(5,2,8),(5,2,10)]
    title_list = ['sim33/36\n POPS PO4 position','DOPX PO4 position','DOPE PO4 position','CHOL ROH position','PPCH PO4 position']
    for sim33_lipid_data_list, sim36_lipid_data_list in zipped_sim33_36_list_of_data:
        sim33_all_frame_numbers_array = numpy.array([array[0] for array in sim33_lipid_data_list[0::3]]) #every 3rd data value starting from first index is a frame number        
        sim36_all_frame_numbers_array = numpy.array([array[0] for array in sim36_lipid_data_list[0::3]]) 
        y_values_array = sim33_lipid_data_list[1] #these are all the same so can just grab one
        sim33_z_values_array = numpy.transpose(numpy.array(sim33_lipid_data_list[2::3])) #a bit of manipulation to get the shape of this array ready for matplotlib
        sim36_z_values_array = numpy.transpose(numpy.array(sim36_lipid_data_list[2::3])) 
        #combine the Z arrays so that sim33 values precede sim36 values
        z_values_array = numpy.concatenate((sim33_z_values_array,sim36_z_values_array),axis=1)
        #combine the frame number arrays after adding in the values for sim36 to avoid overlap with sim33:
        sim36_all_frame_numbers_array += sim33_all_frame_numbers_array[-1]
        x_values_array = numpy.concatenate((sim33_all_frame_numbers_array,sim36_all_frame_numbers_array))/ 10.0
        current_tuple = ax_positions.pop(0)
        ax = fig.add_subplot(current_tuple[0],current_tuple[1],current_tuple[2]) #position plots for each lipid type as you iterate through
        contour_plot = ax.contourf(x_values_array,y_values_array,z_values_array, levels = numpy.logspace(-8.0,6.0,num=50,base=2.0))
        cbar = matplotlib.pyplot.colorbar(contour_plot)
        cbar.ax.set_ylabel('% particles')
        ax.set_title(title_list.pop(0))
#        ax.set_ylabel('Distance from\nlipid centroid (nm)')
        ax.set_ylim((10,45))
        ax.set_xlim((0,1100))
    ax.set_xlabel('Time (ns)')

    fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/leaflet_analysis/stratification_tracking_Dijon.png',dpi=300)
    '''


    #-----------------adjust code below after I get the full pass of data pickled--------------#



#    for lipid_data_list in list_of_data:
#        #have to merge some of the data values together to make amenable to matplotlib contour plotting:
#        all_frame_numbers_array = numpy.array([array[0] for array in lipid_data_list[0::3]]) #every 3rd data value starting from first index is a frame number
#        y_values_array= lipid_data_list[1] #these are all the same so can just grab one
#        z_values_array= numpy.transpose(numpy.array(lipid_data_list[2::3])) #a bit of manipulation to get the shape of this array ready for matplotlib
#        
#        ax = fig.add_subplot(ax_position) #position plots for each lipid type as you iterate through
#        contour_plot = ax.contourf(all_frame_numbers_array,y_values_array,z_values_array, levels = numpy.logspace(-8.0,6.0,num=50,base=2.0)) #plot the data for each timestep individually in the loop
#        ax_position += 1 #place the next lipid species stratification tracking plot below
#        cbar = matplotlib.pyplot.colorbar(contour_plot)
#        cbar.ax.set_ylabel('% particles')
#        ax.set_ylim((10,40))
#        ax.set_title(title_list.pop(0))
#        ax.set_ylabel('Distance from\nlipid centroid (nm)')
#    ax.set_xlabel('Frame number')
#    fig.set_size_inches(8.5,14)
#    fig.savefig('test_stratification_tracking.png',dpi=300)

def positional_probability_virus_species(queue_IPC,universe_object):
    '''Perform the analysis component of the 4D heat map positional probability code for protein species only (as per Mark's request) on a given replicate on a given core.'''
    #set the number of bins early on:
    num_bins = 80 #when I increase to 100+, I have major issues with scene clipping, and I don't know why...
    #as I'm currently using 1000 A ranges in x,y,z direction, this means my bins (cubes) have 12.5 A sides

    #provide a list of selection strings for the molecular species I want to perform positional probability analysis on:
    virion_protein_selection_string_list = ['bynum 1:289680','bynum 289681:338808','bynum 338809:344388'] #all HA, NA, M2 molecules

    #define a function that produces the MDAnalysis selections before trajectory iteration as these are dynamically updated:
    def produce_MDA_selections(list_selection_strings):
        list_selections= []
        for selection in list_selection_strings:
            list_selections.append(universe_object.selectAtoms(selection))
        #now return the ordered list of selections (same order as the strings provided above)
        return list_selections

    selection_list = produce_MDA_selections(virion_protein_selection_string_list)
    #now iterate through the trajectory to obtain the list of numpy arrays (coordinates for specified selections) and histogram these right away to avoid storing large concatenated coordinate arrays [again trying to reduce RAM footprint]       
    frame_counter = 0
    #since Danny's simulation data skips every 10th frame, I'll want to skip every 10 frames on all replicates except for Danny's to maintain balance for comparison
    if universe_object.selectAtoms('resname DUPC').numberOfResidues() > 0:
        frame_skip_value = 1 #don't skip frames as Danny's trajectory is already filtered to every 10th frame
    else: #all other replicates should skip 10 frames for balanced comparison as they are not filtered
        frame_skip_value = 10
    for ts in universe_object.trajectory[::frame_skip_value]: 
        list_coordinate_arrays = [selection.coordinates() for selection in selection_list] #a list of numpy arrays of coordinates in the same order as the selection string lists
        if frame_counter == 0: #initialize list of histogram return objects after the first frame
            histogram_object_list = [numpy.histogramdd(coordinates,bins=(num_bins,num_bins,num_bins),range=[(250,1250),(250,1250),(10,1010)]) for coordinates in list_coordinate_arrays] #will be a list of (histogram, edges) objects with molecular species type ordered as before
            histogram_array = numpy.array([element[0] for element in histogram_object_list]) #values will be added in to this in subsequent frames
            edges_list = [element[1] for element in histogram_object_list] #this won't change in subsequent frames
            frame_counter += 1
        else: #for subsequent frames I'll want to add the species histograms together rather than retaining all the numpy coordinates (again, trying to reduce RAM usage)
            new_histogram_object_list = [numpy.histogramdd(coordinates,bins=(num_bins,num_bins,num_bins),range=[(250,1250),(250,1250),(10,1010)]) for coordinates in list_coordinate_arrays] 
            new_histogram_array = numpy.array([element[0] for element in new_histogram_object_list]) #again, now using numpy array to simplify element-wise addition
            histogram_array += new_histogram_array 
            frame_counter += 1
        print multiprocessing.current_process().name, 'frame counter : ', frame_counter #for tracking multi-core progress at terminal
    #so the cumulative positional histogram for this replicate should be stored in histogram_array

    #I want the midpoints of the bins in x,y,z dimensions for plotting purposes (i.e., the centers of the bins/spheres for heat mapping):
    def bin_edges_to_midpoints(array_bin_edges):
        midpoint_values = []
        for index in range(array_bin_edges.size-1):
            current_midpoint_value = (array_bin_edges.flat[index]+ array_bin_edges.flat[index+1]) /2.0
            midpoint_values.append(current_midpoint_value)
        return numpy.array(midpoint_values)

    #produce a list of shape (num_bins,) 1D numpy arrays of x,y,z midpoints for each molecular species:
    list_midpoints= [] #should have format: [ [x1array,y1array,z1array],[x2array,y2array,z2array],...]
    for bin_edge_array in edges_list:
        x_bin_midpoints,y_bin_midpoints,z_bin_midpoints = [bin_edges_to_midpoints(edge_array) for edge_array in bin_edge_array]
        list_midpoints.append([x_bin_midpoints,y_bin_midpoints,z_bin_midpoints])
    #so should have a nested list of xyz arrays for each of the molecular species with the original ordering preserved

    #return relevant data to the assigned queue object:
    queue_IPC.put([frame_counter,list_midpoints,histogram_array])

def capture_plot_positional_probability_virus_species(Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no'):
    '''Captures the multicore data produced by the similarly named function and plots it for each of the parsed replicates using the PyMOL-based approach. For plotting purposes, I think I'll still have to call this function separately (using pickled data) for each replicate because of the way PyMOL behaves in the scripting context.'''
    if from_pickle == 'no': #if plotting for the first time OR updating the results, pickle the new data to avoid re-parsing for plot modificaitons
        sim35_list_of_data = Q2.get()
        sim36_list_of_data = Q3.get()
        sim37_list_of_data = Q4.get()
        sim38_list_of_data = Q5.get()
        sim39_list_of_data = Q6.get()
        sim_Danny_list_of_data = Q7.get()
        pickle.dump(sim35_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim35_positional_probability_virus_species_March_2014.p','wb'))
        pickle.dump(sim36_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim36_positional_probability_virus_species_March_2014.p','wb'))
        pickle.dump(sim37_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim37_positional_probability_virus_species_March_2014.p','wb'))
        pickle.dump(sim38_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim38_positional_probability_virus_species_March_2014.p','wb'))
        pickle.dump(sim39_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim39_positional_probability_virus_species_March_2014.p','wb'))
        pickle.dump(sim_Danny_list_of_data,open('/sansom/n05/bioc1009/positional_prob_pickle/sim_Danny_positional_probability_virus_species_March_2014.p','wb'))
    elif from_pickle == 'yes': #simply replotting the data from pickled storage
        sim35_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim35_positional_probability_virus_species_March_2014.p','rb'))
        sim36_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim36_positional_probability_virus_species_March_2014.p','rb'))
        sim37_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim37_positional_probability_virus_species_March_2014.p','rb'))
        sim38_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim38_positional_probability_virus_species_March_2014.p','rb'))
        sim39_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim39_positional_probability_virus_species_March_2014.p','rb'))
        sim_Danny_list_of_data = pickle.load(open('/sansom/n05/bioc1009/positional_prob_pickle/sim_Danny_positional_probability_virus_species_March_2014.p','rb'))
    #so, at this stage, regardless of the pickle switch, the appropriate data lists are available for the subsequent PyMOL code
    #I think it will be necessary to run the PyMOL code below separately for each parsed replicate because of the way that program behaves with loop control flow
    def produce_grid_positional_probability_maps(number_bins=80,selected_sim_name='sim35'):
        '''This function should produce a PyMOL-based panel (not strictly using the PyMOL grid command) of positional probability maps for the different protein species. It may also be useful to output multiple viewing angles or cross-sections but I'll start simple for now.'''
        #define some utility functions before I start building up the panel of positional probability maps
        def produce_xyz_file_spherical_bins(xyz_filename,species_x_bin_midpoints,species_y_bin_midpoints,species_z_bin_midpoints):
            '''The purpose of this function is to produce the cubic grid of x,y,z bin centers (cubes) in a .xyz coordinate file. Eventually these cubes may be selectively coloured in a probability heat map in PyMOL. So, I'm basically filling the entire coordinate system around the limits of the virus with small cubes that can be coloured to track the locations of different species of protein.'''
            with open(xyz_filename,'w') as xyz_file:
                #the first line of the xyz file should contain the number of 'atoms:'
                xyz_file.write(str(len(species_x_bin_midpoints) * len(species_y_bin_midpoints) * len(species_z_bin_midpoints)) + '\n')
                #then a comment line:
                xyz_file.write('temporary xyz file\n')
                #produce the cube of cubes by iterating through all possible combinations of bin centers:
                for x in species_x_bin_midpoints:
                    for y in species_y_bin_midpoints:
                        for z in species_z_bin_midpoints:
                            #the first column is the atom symbol; I will use 'C' as filler
                            xyz_file.write('C ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n')

        def produce_pdb_file_spherical_bins(xyz_filename,pdb_filename,histogram_molecular_species):
            '''The purpose of this function is to take an .xyz-formatted file for a given molecular species, which contains a cube of spherical bins, and convert this to a pdb file with b factors set to the histogrammed Cartesian positional probabilities.'''
            print 'starting produce_pdb_file_spherical_bins() function for',sim_name #debugging
            #load the .xyz-formatted file in PyMOL:
            pymol.cmd.load(xyz_filename)
            print 'xyz file loaded'
            #use PyMOL to save the .xyz file as a PDB file:
            pymol.cmd.save(pdb_filename,'all') #save all coords to a new PDB file
            print 'pdb file saved'
            #delete all of the currently loaded PyMOL objects:
            pymol.cmd.delete('all')
            #now load the PDB file that was just produced:
            pymol.cmd.load(pdb_filename)
            print 'pdb file loaded'
            #tricky: put the histogrammed data into a special stored object that is accessible by PyMOL:
            pymol.stored.B_factor_list = histogram_molecular_species.ravel().tolist()
            print 'length B factor list:', len(pymol.stored.B_factor_list)
            #clear out any non-zero B factors in the current PDB:
            pdb_object_name = pdb_filename.replace('.pdb','')
            pymol.cmd.alter(pdb_object_name,"b=0.0")
            #update the B factor values to reflect the positional probabilities at each Cartesian sphere/bin:
            print 'before alter'
            pymol.cmd.alter(pdb_object_name,"b=stored.B_factor_list.pop(0)") #note use of stored object AND iterator requirement
            print 'after alter'
            #now save the updated PDB file with 'B factors' and clear out any objects loaded in PyMOL for good measure:
            pymol.cmd.save(pdb_filename,'all')
            pymol.cmd.delete('all') 

        def place_individual_map_grid(pdb_filename,translation_vector,sphere_vdw_diameter=1000.0/number_bins):
            '''The purpose of this function is to place a given PDB file (with preset B factors for positional probabilities) into the specified panel region (based on translation_vector) in a PyMOL scene, and to setup the representation for rendering. However, I will probably perform the actual final settings and render call in another function.
            sphere_vdw_radius is currently estimated from the outer diameter and the number of bins, but there's probably a more accurate method'''
            #load the PDB of interest (which contains B factors representing positional probabilities):
            pymol.cmd.load(pdb_filename)
            object_name = pdb_filename.replace('.pdb','')
            #hide all of the unnecessary visuals from the current object
            pymol.cmd.hide('everything',object_name)
            #now, all of the empty bins should have b=0, so I only want to display spheres with b > 0:
            pymol.cmd.show('spheres','b > 0') 
            #adjust the VDW radii of the spherical bins:
            pymol.stored.sphere_vdw_radius = sphere_vdw_diameter/2.0
            pymol.cmd.alter(object_name,'vdw=stored.sphere_vdw_radius')
            #since grid mode doesn't seem to be working for me, try to translate objects manually along x,y,z axes:
            pymol.cmd.translate(translation_vector,object_name) #for proper placement in panel
            #now a subsequent function will be used to render all of the panel representations in a single scene
            #I will activate the rainbow spectrum colouring near the end of the workflow, as this seems to work best

        def place_color_bar():
            '''The purpose of this function is to produce and position a color bar in the PyMOL figure to show the value scale of the heat map. Note that I'm currently using a scale from 0.0 to 1.0 but my values actually go up to 10 (this seems to improve contrast for smaller values).'''
            #produce an xyz file with about 100 coordinate values, each separated horizontally by ~20A
            x_values = range(0,2000,20)
            #use constant y and z values for a straight line color bar; may need some adjusting to get the position right;
            y_value = 500; z_value = 0 
            #if sim_name == 'sim33':
                #y_value = -350 #shift cb down a bit for the case of the larger equilibrating vesicles
            with open('color_bar.xyz','w') as xyz_file:
                #the first line should cotain the total number of atoms:
                xyz_file.write(str(len(x_values)) +'\n')
                #the second line is to be used for commenting:
                xyz_file.write('colorbar .xyz file\n')
                for x_value in x_values: #write in all of the coords, including filling the first column with 'Carbon atoms'
                    xyz_file.write('C ' + str(x_value) + ' ' + str(y_value) + ' ' + str(z_value) + '\n')
            #load the xyz file in PyMOL:
            pymol.cmd.load('color_bar.xyz') 
            #save the xyz file as a PDB file and delete the xyz object:
            pymol.cmd.save('color_bar.pdb','color_bar')
            pymol.cmd.delete('color_bar')
            #load the PDB file and display it as a color bar:
            pymol.cmd.load('color_bar.pdb')
            pymol.cmd.hide('everything','color_bar')
            #the B values will range over 0 to 1 with an even spacing over the length of the number of 'coords:'
            pymol.stored.color_bar_B_values = numpy.linspace(0.0,200.0,num=len(x_values)).ravel().tolist()
            #use the progressively increasing list of artificial B values to produce the gradations in the colorbar colors by subtituting the B values and using the rainbow scheme:
            pymol.cmd.alter('color_bar',"b=0.0")
            pymol.cmd.alter('color_bar',"b=stored.color_bar_B_values.pop(0)") #note use of stored object AND iterator requirement
            pymol.cmd.save('color_bar.pdb', 'color_bar') 
            #alright delete the color_bar object and load the PDB file with correct B values, as it seems to save correctly:
            pymol.cmd.delete('color_bar')
            pymol.cmd.load('color_bar.pdb')
            pymol.cmd.alter('color_bar','vdw=32.0')
            pymol.cmd.show('spheres','color_bar') #this line seems to cause a segfault when rendering if I zoom out too far?!
            #hold off on activating the rainbow/color scheme until later in the workflow as this seems to work best

        def place_decorations():
            '''Should place labels, etc., where they are needed in the PyMOL scene. This will probably make extensive use of pseudoatom objects. Too many problems getting labels to render with povray so will shelf this function and put labels in post-production.'''
            def label_species_heat_maps(label_name,object_name):
                '''Use pseudoatom objects to place labels above the virion heat maps for each lipid species.'''
                pseudoatom_object_name = object_name.replace('temp','pseudoatom')
                #create a pseudoatom object anchored at the center of the specified object:
                pymol.cmd.pseudoatom(pseudoatom_object_name,object_name)
                #now translate the pseudoatom to a position above the virion heat map for the given species:
                pymol.cmd.translate([0,350,300],pseudoatom_object_name)
                #hide the pseudoatom particle, as I don't need that:
                pymol.cmd.hide('everything',pseudoatom_object_name)
                #place the label at the pseudoatom:
                #pymol.cmd.label(pseudoatom_object_name,label_name)
                pymol.cmd.label(pseudoatom_object_name,'"test_label"')
                pymol.cmd.show('spheres',pseudoatom_object_name)
                pymol.cmd.show('label',pseudoatom_object_name)
                pymol.cmd.alter(pseudoatom_object_name,'vdw=16.0')
            #now call the above function to place species labels above each of the heat maps:
            for object_name in ['all_HA','all_NEU','all_M2']:
                label_species_heat_maps('"%s"' % object_name.replace('_temp',''),object_name) #note the strict string quoting requirements for PyMOL
    
        def render_scene(simulation_name,final_frame):
            '''The purpose of this function is to adjust final render settings and produce the output image file of the positional probability maps on a panel. For now, I can't get a surface representation to work, but this is among the things I'd like to improve for this code.'''
            #for sims#35/36 I'll want to rotate the HA1 + PPCH PO4 representation so that the Z axis points up along Y, so that it is easier to see that HA1 is 'shooting off' the sphere in the correct orientation
            #I'll want a white background:
            pymol.cmd.bg_color('white')
            #put the spectrum colouring command near the end and use explicit limits to ensure that everything is on the same scale
            #for now I"m adjusting the max value below to scale based on manual inspection of the max values in the replicate histograms as there is currently no normalization to 1 above
            pymol.cmd.spectrum('b','rainbow',minimum=0.0,maximum=200) #testing a pseudo-non-linear scale that goes up to 10, but colouring stops changing after 1 (red)
            #the following surface code most likely is not working because my particles (spherical bins) are too far apart, and trying to produce smaller spheres (by increasin num_bins) causes massive scene clipping/cropping for some reason
            #CAREFUL: my PDB spherical bins are all based on HETATM entries which are, by default, hidden when loaded, so I need to remove the hiding effect in order to draw a surface:
            #pymol.cmd.flag('ignore','all','clear')
            #pymol.cmd.rebuild()
            #try employing a surface representation:
            #pymol.cmd.create('special_obj', 'POPS_temp and b > 0') 
            #pymol.cmd.ramp_new('ramp_obj', 'special_obj', [0, 1], 'rainbow')
            #pymol.cmd.set('surface_color', 'ramp_obj', 'special_obj')
            #pymol.cmd.show('surface','special_obj')
            pymol.cmd.rebuild()
            #I've used this setting for nice renders in the past:
            pymol.cmd.set('antialias','2')
            #try to avoid ray tracing segmentation faults by controlling the amount of memory that PyMOL uses for ray tracing
            # memory usage depends on hash_max^3; hash_max = 240 --> ~243 MB hash + data
            pymol.cmd.set('hash_max',200)
            #adjust zoom appropriately:
            #if sim_name == 'sim33':
                #pymol.cmd.zoom('POPS_temp or DOPX_temp or DOPE_temp or PPCH_temp or CHOL_temp or Z2655_temp',100) #zoom excludes protein selections when working with vesicle sim
            #elif sim_name == 'sim38':
                #pymol.cmd.zoom('POPS_temp or DOPX_temp or DOPE_temp or PPCH_temp or CHOL_temp or Z2655_temp or first_HA or all_HA or all_NEU or all_M2 or FORS_temp',100)
            #else: #but zoom can include protein selections for virion replicates
            pymol.cmd.zoom('all_HA or all_NEU or all_M2',100)

            #render the image and quit PyMOL:
            print 'before ray'
            pymol.cmd.ray(width=2000,height=2000,renderer=1) #render with povray rather than internal renderer; deal with DPI stuff after as png command not doing what I want
            #the above external render will produce a file named tmp_pymol.png, which isn't ideal, but at least it's doing a nice ray trace without segfaulting, which was a problem with internal ray tracing
            #also doesn't include text labels, but I've had issues with the png command after using the external renderer so I'll just have to deal with labels in Photoshop
            print 'after ray'
            #okay, I'll want to rename the rendered file (which is always named 'tmp_pymol.png') to something informative
            #sim_name = 'sim35' sim_name will actually be specified upon calling in multicore version
            new_png_filename = simulation_name + '_' + 'frame_' + str(final_frame) + 'positional_probability.png'
            os.rename('tmp_pymol.png',new_png_filename)
            pymol.cmd.quit()
        #---------------------------------------------------------------------------------------------#            
        #end of utility function definitions/ start of main function work
        
        #to reduce physical memory footprint, iterate through the replicates and load/ plot them one at a time:
        for sim_name,list_of_data in zip(['sim35','sim36','sim37','sim38','sim39','sim_Danny'],[sim35_list_of_data,sim36_list_of_data,sim37_list_of_data,sim38_list_of_data,sim39_list_of_data,sim_Danny_list_of_data]):
            if not sim_name == selected_sim_name: continue #seem to be limited to plotting one heat map panel at a time
            print 'starting plot process for', sim_name
            #perform some necessary commands to initiate PyMOL:
            import __main__
            __main__.pymol_argv = ['pymol','-c'] # Pymol: quiet and no GUI (currently just no GUI as the -q flag masks any errors)
            import pymol
            pymol.finish_launching()
            #step 1: produce the xyz cube grid file for each of the molecular species I'm going to map-------------------------
            xyz_files = ['all_HA.xyz','all_NEU.xyz','all_M2.xyz']
            #if sim_name == 'sim33':
                #xyz_files = xyz_files[:-4] #exclude protein representations for vesicle sim
            #if sim_name == 'sim38': 
                #xyz_files = xyz_files + ['FORS_temp.xyz']
            for xyz_filestring,centers in zip(xyz_files,list_of_data[1]):
                #produce all the xyz files for each of the species of interest:
                produce_xyz_file_spherical_bins(xyz_filestring,centers[0],centers[1],centers[2])
            print 'xyz files produced for', sim_name
            #end step 1----------------------------------------------------------------------------------------
            #step 2: produce the pdb files with appropriate B factor columns---------------------------------------------------
            for xyz_file,species_histogram in zip(xyz_files,list_of_data[2]):
                produce_pdb_file_spherical_bins(xyz_file,xyz_file.replace('xyz','pdb'),species_histogram)
            print 'pdb files produced for', sim_name
            #end step 2--------------------------------------------------------------------------------------------------------
            #step 3: place each PDB representation at the appropriate grid location and prepare for rendering-----------------
            pdb_filenames = [xyz_filename.replace('xyz','pdb') for xyz_filename in xyz_files]
            translation_vectors = [[-700,400,0],[300,400,0],[1300,400,0]] #translate heat maps so they can be seen individually in the panel
            #if sim_name == 'sim33': #need different spacing for large vesicle equilibration stages
                #translation_vectors = [[0,400,0],[900,400,0],[1800,400,0],[0,-500,0],[900,-500,0],[1800,-500,0]]
            #elif sim_name == 'sim38': #there's an extra plot for FORS here:
                #translation_vectors.append([3000,-500,0])
            for pdb_file in pdb_filenames:        
                place_individual_map_grid(pdb_file,translation_vectors.pop(0),sphere_vdw_diameter=1000.0/number_bins)
            print 'plot elements placed on panel for',sim_name
            #end step 3----------------------------------------------------------------------------------------------------
            #step 4: place an appropriate color bar (unfortunately I can't get labels to work well with external renderer; png command causes issues)
            print 'before color bar'
            place_color_bar()
            print 'after color bar'
            #end step 4-------------------------------------------------------------------------------------------
            #step 5: finally, render the grid of heat maps-----------------------------------------------------------------
            print 'before render call'
            render_scene(sim_name,list_of_data[0])
            print 'after render call' #control flow stalling here without error?!
            #end step 5--------------------------------------------------------------------------------------------------
            
        
    #call the above function to get the heat maps rendered:
    #the problem I have right now is that PyMOL won't 'restart' on each iteration of the rendering loop in this plotting code:
    #so, basically, only the selected element gets plotted
    produce_grid_positional_probability_maps(selected_sim_name = 'sim39')
    print 'exiting plot/capture function... (shoud finish)'


def protein_proximity_analysis(queue_IPC,universe_object):
    '''Multicore-adapted MDA code to parse the % of a given lipid species around a given protein species in a given virion simulation condition. Currently under development as of December 31/ 2012. This code will not distinguish between leaflets, at least for now. It would probably also be sensible to normalize the data for the protein populations (because there is much more HA than M2, etc.), so I will probably divide the count or percent arrays by the total number of proteins of a given type so that I get % lipid population/ protein on average.'''
    #adjusting code for Biophysics 2013 poster requirements:
    #sims 35,36,and 37 should all have the same topology/index numbers for protein selections, so each core can use consistent code
    #sim 38 should also be the same for protein selections, but obviously has the new FORS glycolipid as well

    #it would be a waste of time to use MDA to grab the total residue counts as I already know these from the topology, so hard code a dictionary for this:
    lipid_residue_count_dictionary = {'POPS':6288, 'DOPX':3753, 'DOPE':2111, 'CHOL':22829, 'PPCH':7931}
    #however, the above dictionary would only really be needed to calculate % values, which I'm less interested in for Biophysics, AND it would be a slightly different dictionary for sim38 with ~68% of PPCH replace with FORS

    #define a generic function that will calculate the % of a given lipid species around each of the three protein types in a given time step:
    #for Biophysics analysis, I simply want this to do raw counts (normalized for protein, but NOT lipid populations)
    def lipid_protein_proximity_calculator(lipid_resname,total_residues_this_type,proximity_threshold):
        '''proximity_threshold in A units for lipid around protein inclusion distance; total_residues_this_type should be an integer representing the total number of lipid_resname species in the system, and will be used for calculating percent values; '''
        HA_based_selection_string = 'resname {lipid_name} and (around {threshold} bynum 1:289680)'.format(lipid_name = lipid_resname, threshold = str(proximity_threshold))
        NA_based_selection_string = 'resname {lipid_name} and (around {threshold} bynum 289681:338808)'.format(lipid_name = lipid_resname, threshold = str(proximity_threshold))
        M2_based_selection_string = 'resname {lipid_name} and (around {threshold} bynum 338809:344388)'.format(lipid_name = lipid_resname, threshold = str(proximity_threshold))
        #now, produce MDA selection objects for each of the 3 protein-lipid proximity tests
        HA_proximal_selection, NA_proximal_selection, M2_proximal_selection = [universe_object.selectAtoms(proximity_string) for proximity_string in [HA_based_selection_string,NA_based_selection_string,M2_based_selection_string]]
        #the number of residues in each of the above AtomGroups serves as a count of the number of lipids of this type within the distance cutoff for a given protein_type:
        species_count_near_HA, species_count_near_NA,species_count_near_M2 = [atom_group.numberOfResidues() for atom_group in [HA_proximal_selection, NA_proximal_selection, M2_proximal_selection]]
        #convert the protein-proximal lipid species counts to % of the total number of residues of this species:
        #species_percent_near_HA, species_percent_near_NA, species_percent_near_M2 = map(lambda x: x/float(total_residues_this_type) * 100.0, [species_count_near_HA, species_count_near_NA,species_count_near_M2]) 
        #return a tuple of the % lipid species near protein values:
        #actually, for Biophysics 2013 analysis, simply return a tuple of the raw count near protein values:
        return (species_count_near_HA, species_count_near_NA, species_count_near_M2)

    #create another dictionary object with the same key names for appending the % proximal counts as I iterate through a trajectory:
    #percent_species_count_dictionary = {'POPS':[],'DOPX':[],'DOPE':[],'CHOL':[],'PPCH':[]}
    raw_species_count_dictionary = {'POPS':[],'DOPX':[],'DOPE':[],'CHOL':[],'PPCH':[]}

    for ts in universe_object.trajectory[::100]: #every N frames for now
        #iterate through all the different lipid species and do the protein-proximal counting/appending of values
        for residue_name in raw_species_count_dictionary.iterkeys():
            species_count_near_HA, species_count_near_NA, species_count_near_M2 = lipid_protein_proximity_calculator(residue_name,lipid_residue_count_dictionary[residue_name],proximity_threshold=12.0) #currently using a 12.0 A cutoff
            #print multiprocessing.current_process().name, 'debugging, species_percent_near_HA:', species_percent_near_HA
            raw_species_count_dictionary[residue_name].append([species_count_near_HA, species_count_near_NA, species_count_near_M2]) #so each of the lists in count_species_count_dictionary can eventually be converted to numpy arrays with shape (N_frames, 3) where the first column of data is for HA, the second column of data is for NA, etc.
            #I will probably divide the respective values by the total number of proteins (to normalize for relative protein abundance in virion), but this can be done in the capture/plot function after the multicore component of the code completes
        #when all of the counting and appendng has been done for the current frame on a given processor, record the progress in the terminal:
        print multiprocessing.current_process().name, 'frame : ', ts.frame
    #just use counts, not %s for Biophysics 2013 analysis:
    queue_IPC.put(raw_species_count_dictionary) #return the dictionary object with the nested percentage lists to the queue

def capture_plot_protein_proximity_analysis(Q2,Q3,Q4,Q5,from_pickle = 'no'):
    '''Capture, modify, and plot the data from the matching multicore analysis function above. As the matching function only works sensibly in the presence of viral proteins, at present sims 35,36 and 37 apply, bit sim 33 does not (so use queue objects 2,3,4). Now adjusting to include queue object 5 as sim38 is available as well.'''
    #if not plotting from pickle files then generate objects directly from the multiprocessing queues and then pickle for safety/future use:
    if from_pickle == 'no':
        #capture the dictionary objects sent back to queues from each core
        sim_35_raw_species_count_dictionary= Q2.get()
        sim_36_raw_species_count_dictionary= Q3.get()
        sim_37_raw_species_count_dictionary= Q4.get()
        sim_38_raw_species_count_dictionary= Q5.get()
        #pickle the dictionaries to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim_35_raw_species_count_dictionary,open('sim_35_raw_species_count_dictionary.p','wb'))
        pickle.dump(sim_36_raw_species_count_dictionary,open('sim_36_raw_species_count_dictionary.p','wb'))
        pickle.dump(sim_37_raw_species_count_dictionary,open('sim_37_raw_species_count_dictionary.p','wb'))
        pickle.dump(sim_38_raw_species_count_dictionary,open('sim_38_raw_species_count_dictionary.p','wb'))
        print 'pickling complete' #debugging
    elif from_pickle == 'yes': #for this case I simply retrieve the dictionaries from the pickled files
        list_dictionaries = [pickle.load(open(x,'rb')) for x in ['sim_35_raw_species_count_dictionary.p','sim_36_raw_species_count_dictionary.p','sim_37_raw_species_count_dictionary.p','sim_38_raw_species_count_dictionary.p']]
        sim_35_raw_species_count_dictionary,sim_36_raw_species_count_dictionary,sim_37_raw_species_count_dictionary,sim_38_raw_species_count_dictionary = list_dictionaries
    #plotting code here
    import matplotlib.pyplot
    matplotlib.rc('xtick', labelsize=18) 
    matplotlib.rc('ytick', labelsize=18) 
    #write a generic plotting function that I can call for each of the simulation conditions:
    def plot_protein_proximity_data(raw_species_count_dictionary,plot_title,outfile_name):
        fig = matplotlib.pyplot.figure()
        ax_HA = fig.add_subplot(311) #HA data in first plot row
        ax_NA = fig.add_subplot(312) #NA data in second plot row
        ax_M2 = fig.add_subplot(313) #M2 data in third plot row
        for residue_name in raw_species_count_dictionary.iterkeys():
            array_of_raw_species_counts = numpy.array(raw_species_count_dictionary[residue_name]) #as described above, HA, NA, and M2 data are stored in columns 1, 2, and 3 respectively
            #plot the data for this lipid residue in each of the 3 subplots (one for each flu protein type):
            protein_num = 0
            for protein_type, normalization_factor in [[ax_HA,80],[ ax_NA,12], [ax_M2,15]]: #divide count arrays by normalization factor, which is simply the abundance of a protein in the virion (otherwise it's not a fair comparison between the different protein types as HA is far more abundant)
                protein_type.plot(numpy.arange(0,array_of_raw_species_counts[...,protein_num].size)/100.0,array_of_raw_species_counts[...,protein_num]/normalization_factor,linewidth=4.0,label=residue_name) 
                protein_num += 1 #increment the column of data used, which corresponds to data for the different protein types
        ax_M2.set_xlabel('Time ' + r'($\mu$s)',fontsize=18) #place x-axis label under bottom plot only
        #ax_HA.set_title(plot_title) #place plot title on top of top row plot
        #ax_HA.set_ylabel('raw lipid species\n within 12 ' + r'$\AA$' + ' of HA \n(population normalized)')
        #ax_NA.set_ylabel('raw lipid species\n within 12 ' + r'$\AA$' + ' of NA \n(population normalized)')
        #ax_M2.set_ylabel('raw lipid species\n within 12 ' + r'$\AA$' + ' of M2 \n(population normalized)')
        ax_HA.set_ylabel('lipid count \nnear HA',fontsize=16)
        ax_NA.set_ylabel('lipid count \nnear NA',fontsize=16)
        ax_M2.set_ylabel('lipid count \nnear M2',fontsize=16)
        for ax in [ax_HA,ax_NA,ax_M2]: ax.set_ylim((0.0,70.0)) #use same y-scale for the 3 protein species plots
        for ax in [ax_HA,ax_NA,ax_M2]: ax.set_xlim((0.0,5.0)) #use same x-scale for the 3 protein species plots, and between sim conditions
        #place a legend below the bottom (M2) axis:
        ax_M2.legend(bbox_to_anchor=(1.17, 3.4), borderaxespad=0.)
        fig.set_size_inches(11.5,6.8)
        matplotlib.pyplot.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)
        fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/protein_percent_species_counts/' + outfile_name, dpi = 300)

    #produce the plots for each of the simulation conditions using the above generic plotting function:
    for raw_species_count_dictionary,plot_title,outfile_name in [(sim_35_raw_species_count_dictionary,'Influenza virion CG-MD simulation at 295K with proteins unrestrained\n (sim35;  every 10 ns parsed; 10 hours parallel)\n','sim35_protein_proximity_analysis_SBCB_seminar.png'),(sim_36_raw_species_count_dictionary,'Influenza virion CG-MD simulation at 295K with proteins restrained\n (sim36;  every 10 ns parsed; 10 hours parallel)\n','sim36_protein_proximity_analysis_SBCB_seminar.png'),(sim_37_raw_species_count_dictionary,'Influenza virion CG-MD simulation at 323K with proteins unrestrained\n (sim37;  every 10 ns parsed; 10 hours parallel)\n','sim37_protein_proximity_analysis_SBCB_seminar.png'),(sim_38_raw_species_count_dictionary,'sim38','sim38_protein_proximity_analysis_SBCB_seminar.png')]:
        plot_protein_proximity_data(raw_species_count_dictionary,plot_title,outfile_name)

def protein_proximity_analysis_v2(queue_IPC,universe_object):
    '''March 2013: developing a new code for analysing the proximity of lipid species to proteins by counting the number of proteins of a given type (HA, NA or M2) that are within a certain cutoff of each molecule of the lipid type in question. I'm going to start off my testing on cholesterol, but may work up to assessing the full set of lipids in the various influenza virions. I think it's important to produce plots at various intervals during a replicate simulation as I want to know if the protein-lipid interactions are CHANGING--otherwise, there isn't much evidence that adding proteins to the vesicle system has had a substantial impact on lipid domains. For now, however, I think it can still be interesting to parse entire trajectories and compare the different replicates to see if there are overall differences.'''
              

    #---------------------------------------------------------------------------------------
    #write a more general function for counting the number of proximal proteins of a given type to a set of lipid residues of a given type:
    def lipid_residue_protein_proximity_counter(lipid_type, protein_type, distance_cutoff = 12.0):
        '''For each of the residues for a given input lipid type and protein type, assess the number of proteins of that type that fall within a specified cutoff distance, and return the set of all these distances in a numpy array. Note that my pre-filter strategy currently does not count lipids that are not within the distance_cutoff of any of the specified proteins--so zeros are not counted for binning at the moment, which may or may not be something I want to change, but I'm trying to cut down on the amount of internal for looping for speed.
        protein_type can be: 'HA', 'NA', or 'M2'
        distance_cutoff is specified as a float in Angstrom units (defaults to 12.0)
        lipid_type can be 'POPS', 'DOPX', 'DOPE', 'CHOL', 'PPCH', 'FORS'
        '''
        #the number of CG particles in the protein is dictated by protein_type, as is the number of proteins
        #now modifying code to use numpy arrays more extensively and avoid MDA looping selections/ distance checks
        if protein_type == 'HA':
            protein_particles = 3621
            num_proteins = 80
            initial_topology_index = 1
            #perform a single MDA selection to pull out all of the HA coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 1:289680').coordinates()
            #print multiprocessing.current_process().name,'shape of protein_coord_array:', protein_coord_array.shape, '(should be: 289680,3)'
            pre_filter_string = 'and (around {threshold} bynum 1:289680)'.format(threshold = str(distance_cutoff))
        elif protein_type == 'NA':
            protein_particles = 4094
            num_proteins = 12
            initial_topology_index = 289681
            #perform a single MDA selection to pull out all of the NA coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 289681:338808').coordinates()
            #print multiprocessing.current_process().name,'shape of protein_coord_array:', protein_coord_array.shape, '(should be: 49128,3)'
            pre_filter_string = 'and (around {threshold} bynum 289681:338808)'.format(threshold = str(distance_cutoff))
        elif protein_type == 'M2':
            protein_particles = 372
            num_proteins = 15
            initial_topology_index = 338809
            #perform a single MDA selection to pull out all of the M2 coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 338809:344388').coordinates()
            #print multiprocessing.current_process().name,'shape of protein_coord_array:', protein_coord_array.shape, '(should be: 5580,3)'
            pre_filter_string = 'and (around {threshold} bynum 338809:344388)'.format(threshold = str(distance_cutoff))


        #as was done with the proteins, I want to obtain a numpy array of coordinates for the all the lipids of the specified lipid type
        specified_lipid_selection_string = 'resname {lipid_residue_name}'.format(lipid_residue_name = lipid_type) 
        lipid_selection = universe_object.selectAtoms(specified_lipid_selection_string)
        lipid_coord_array = lipid_selection.coordinates()
        #print multiprocessing.current_process().name, 'shape of lipid_coord_array:', lipid_coord_array.shape, '(should be: number of lipid coords: 3)'
        #I'm also going to want to know the number of particles than are within a single lipid residue for striding purposes:
        num_lipids = float(len(lipid_selection.resids()))
        particles_per_lipid = float(lipid_coord_array.shape[0]) / num_lipids  
        #print multiprocessing.current_process().name, 'number of particles per lipid:', particles_per_lipid

        #now, I want a function that can 'stride' through the protein and lipid numpy coordinate arrays to average coordinate sets into the centers of geometry for individual proteins or lipid residues into new numpy arrays
        def produce_COG_arrays(particles_per_molecule, coordinate_array, num_molecules, molecule_type):
            '''Given a numpy coordinate_array of shape (num_particles, 3), use the particles_per_molecule to slice through the numpy array, average the Cartesian coordinates for each molecule/assembly to produce the COG of that assembly, and return an ordered numpy array of the COGs with shape (num_molecules, 3).
            molecule_type can be 'HA','NA','M2', or 'lipid'; the point is to only use TMDs for the protein COG calculations.'''
            COG_list = [] #append molecule/unit COG coords here, in topological order
            remaining_molecules = num_molecules 
            start_row_coefficient = 0
            def TMD_particle_selector(input_array):
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

            while remaining_molecules > 0:
                start_row = start_row_coefficient * particles_per_molecule
                end_row = ((start_row_coefficient + 1) * particles_per_molecule) 
                current_molecule_coordinate_array = coordinate_array[start_row:end_row,...]
                #deal with TMD selection if applicate (lipid cases will simply return the same array)
                current_molecule_coordinate_array = TMD_particle_selector(current_molecule_coordinate_array)
                #assert current_molecule_coordinate_array.shape == (particles_per_molecule,3), "current_molecule_coordinate_array has incorrect shape"
                #average together the above set of Cartesian coordinates for this molecule to obtain its COG:
                current_molecule_COG_array = numpy.average(current_molecule_coordinate_array,axis = 0)
                assert current_molecule_COG_array.shape == (3,), "current_molecule_COG_array has incorrect shape"
                #append the COG for the current molecule to the list and then increment/decrement varibles and continue loop:
                COG_list.append(current_molecule_COG_array)
                start_row_coefficient += 1 #pattern for row limits is 0:n , n:2n , 2n:3n , etc.
                remaining_molecules -= 1
            #convert the COG_list to a numpy array and verify its shape:
            COG_array = numpy.array(COG_list)
            assert COG_array.shape == (num_molecules, 3), "COG_array has incorrect shape."
            #finally return the numpy array containing the ordered set of per-molecule or per-assembly center of geometry values:
            return COG_array
    
        #use the above function to obtain two ordered per-molecule or per-assembly COG coordinate arrays:
        lipid_COG_coordinate_array = produce_COG_arrays(particles_per_lipid, lipid_coord_array, num_lipids, 'lipid')
        protein_COG_coordinate_array = produce_COG_arrays(protein_particles, protein_coord_array, num_proteins, protein_type)

        #I want to calculate the distance between each lipid of the specified type and each protein of the specified type
        import scipy, scipy.spatial.distance
        distance_matrix = scipy.spatial.distance.cdist(lipid_COG_coordinate_array,protein_COG_coordinate_array,'euclidean')
        #save distance_matrix to a text file for debugging purposes:
        numpy.savetxt('distance_matrix_debug.txt',distance_matrix,newline ="\n")
        #another potentially useful debugging value is the closest lipid-protein approach in the distance_matrix:
        print 'closest lipid-protein approach (A):',distance_matrix.min()
        assert distance_matrix.shape == (lipid_COG_coordinate_array.shape[0],protein_COG_coordinate_array.shape[0]), "distance_matrix has incorrect shape."
        #so, each row of the distance_matrix should correspond to a lipid and each column in that row should correspond to the distance between that lipid and one of the proteins of the specified type

        #now I want to use numpy to count the number of occurrences of lipid-protein COG-COG contacts within the distance_matrix that fall within the distance cutoff, and retain the per-lipid (per-row) total counts of proteins within the cutoff
        array_per_lipid_protein_proximal_counts = (distance_matrix < distance_cutoff).sum(axis=1) #sum across numpy Boolean arrays (1 = True [within cutoff]; 0 = False)
        assert array_per_lipid_protein_proximal_counts.shape == (lipid_COG_coordinate_array.shape[0],), "array_per_lipid_protein_proximal_counts has incorrect shape."
        #return the array containing the ordered set of per-lipid residue protein-proximal counts
        return array_per_lipid_protein_proximal_counts
    #---------------------------------------------------------------------------------------

    range_list = [[0,10],[0,10]] #the xmin,xmax and ymin,ymax protein proximity count values for binning (so currently using a cap of 10 proteins within the cutoff as this seems like a difficult number to exceed)
    #note that, at the moment, the code doesn't iterate through a trajectory, it will simply work on the first frame (this is intentional for testing purposes)
    #now start modifying code to iterate over a trajectory and produce histogrammed results (but still working with a small test data set)
    list_count_arrays = []
    for ts in universe_object.trajectory[::100]: #parse every Nth frame
        if ts.frame == 1: #for the first frame only
            for protein_type in ['HA','NA','M2']:
                list_count_arrays.append(lipid_residue_protein_proximity_counter('CHOL',protein_type,distance_cutoff=120.0))
            #now, I want to generate a 2D histogram for HA relative to NA, and another 2D histogram for HA relative to M2:
            HA_NA_histogram, HA_edges_1, NA_edges  = numpy.histogram2d(list_count_arrays[0],list_count_arrays[1],bins=11,range=range_list)
            HA_M2_histogram, HA_edges_2, M2_edges = numpy.histogram2d(list_count_arrays[0],list_count_arrays[2],bins=11,range=range_list)
            print multiprocessing.current_process().name, 'frame: ', ts.frame
        else: #for all subsequent frames
            list_count_arrays = [] #reset for reuse 
            for protein_type in ['HA','NA','M2']:
                list_count_arrays.append(lipid_residue_protein_proximity_counter('CHOL',protein_type,distance_cutoff=120.0))
            HA_NA_histogram_current, HA_edges_current, NA_edges_current = numpy.histogram2d(list_count_arrays[0],list_count_arrays[1],bins=11,range=range_list)
            HA_M2_histogram_current, HA_edges_2_current, M2_edges_current = numpy.histogram2d(list_count_arrays[0],list_count_arrays[2],bins=11,range=range_list)
            #I think it would not be good if the bins / edge list values changed during the trajectory (as summing histogram frequencies over the length of a trajectory would then become invalid), so check for this:
            assert numpy.subtract(HA_edges_current,HA_edges_1).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(NA_edges_current,NA_edges).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(HA_edges_2_current,HA_edges_2).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(M2_edges_current,M2_edges).sum() == 0, "histogram edge inconsistency"
            #if the edge list assertions pass, then combine histograms:
            HA_NA_histogram = numpy.add(HA_NA_histogram,HA_NA_histogram_current)
            HA_M2_histogram = numpy.add(HA_M2_histogram,HA_M2_histogram_current)
            print multiprocessing.current_process().name, 'frame: ', ts.frame

    queue_IPC.put([[HA_NA_histogram,HA_edges_current, NA_edges_current],[HA_M2_histogram,HA_edges_2_current, M2_edges_current]]) #return the object to the IPC queue




def capture_plot_protein_proximity_analysis_v2(Q2,Q3,Q4,Q5,from_pickle = 'no'):
    '''Capture, modify, and plot the data from the matching multicore analysis function above. As the matching function only works sensibly in the presence of viral proteins, at present sims 35,36,37, and 38 apply, but sim 33 does not (so use queue objects 2,3,4,5).'''
    #if not plotting from pickle files then generate objects directly from the multiprocessing queues and then pickle for safety/future use:
    if from_pickle == 'no':
        #capture the lists of 2d histograms from the IPC queue objects, with format: [[HA_NA_histogram,HA_edges_current, NA_edges_current],[HA_M2_histogram,HA_edges_2_current, M2_edges_current]]
        sim_35_list_2d_histograms= Q2.get()
        sim_36_list_2d_histograms= Q3.get()
        sim_37_list_2d_histograms= Q4.get()
        sim_38_list_2d_histograms= Q5.get()
        #pickle the lists to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim_35_list_2d_histograms,open('sim_35_list_2d_histograms.p','wb'))
        pickle.dump(sim_36_list_2d_histograms,open('sim_36_list_2d_histograms.p','wb'))
        pickle.dump(sim_37_list_2d_histograms,open('sim_37_list_2d_histograms.p','wb'))
        pickle.dump(sim_38_list_2d_histograms,open('sim_38_list_2d_histograms.p','wb'))
        print 'pickling complete' #debugging
    elif from_pickle == 'yes': #for this case I simply retrieve the lists from the pickled files
        list_of_lists= [pickle.load(open(x,'rb')) for x in ['sim_35_list_2d_histograms.p','sim_36_list_2d_histograms.p','sim_37_list_2d_histograms.p','sim_38_list_2d_histograms.p']]
        sim_35_list_2d_histograms, sim_36_list_2d_histograms, sim_37_list_2d_histograms, sim_38_list_2d_histograms = list_of_lists
    #with the histogram data now available, whether by direct retrieval from queue objects or by retrieval from pickled objects, go ahead and write some test plotting code to visualize the data:
    #write a general function that can be called to plot the data for each replicate:
    import matplotlib.pyplot
    def replicate_contour_plotter(lipid_name,replicate_name,histogram_object_list,outfile_name,cutoff):
        '''replicate_name can be 'sim35','sim36','sim37', or 'sim38' 
        lipid_type can be 'POPS', 'DOPX', 'DOPE', 'CHOL', 'PPCH', 'FORS' '''
        fig = matplotlib.pyplot.figure()
        #testing: plot lipid species protein-proximity data for HA/NA and HA/M2 side-by-side
        ax_HA_NA = fig.add_subplot(121) #HA_NA data in first plot column 
        ax_HA_NA.set_xlabel('NA proteins within {cutoff_string} $\AA$'.format(cutoff_string = cutoff))
        ax_HA_NA.set_ylabel('HA proteins within {cutoff_string} $\AA$'.format(cutoff_string = cutoff))
        #matplotlib.pyplot.title(50 * ' ' + '{lipid_string} protein proximity counting ({replicate_string} test case)'.format(replicate_string = replicate_name, lipid_string = lipid_name))
        ax_HA_M2 = fig.add_subplot(122) #HA_M2 data in second plot column 
        ax_HA_M2.set_xlabel('M2 proteins within {cutoff_string} $\AA$'.format(cutoff_string = cutoff))
        ax_HA_M2.set_ylabel('HA proteins within {cutoff_string} $\AA$'.format(cutoff_string = cutoff))
        histogram_list = [histogram_object_list[0][0],histogram_object_list[1][0]]
        #obtain the extents from the queue/pickled objects as well:
        xedges_HA_NA = histogram_object_list[0][1]
        #print 'xedges_HA_NA:',xedges_HA_NA #debug print to see what's going on
        yedges_HA_NA = histogram_object_list[0][2]
        xedges_HA_M2 = histogram_object_list[1][1]
        yedges_HA_M2 = histogram_object_list[1][2]
        extent_list = [[yedges_HA_NA[0],yedges_HA_NA[-1],xedges_HA_NA[0],xedges_HA_NA[-1]],[yedges_HA_M2[0],yedges_HA_M2[-1],xedges_HA_M2[0],xedges_HA_M2[-1]]]
        contour_plot_1 = ax_HA_NA.contourf(histogram_list[0],extent=extent_list[0],origin='lower',levels = numpy.linspace(0,numpy.amax(histogram_list[0])))
        contour_plot_2 = ax_HA_M2.contourf(histogram_list[1],extent=extent_list[1],origin='lower',levels = numpy.linspace(0,numpy.amax(histogram_list[1])))
        cbar_1 = matplotlib.pyplot.colorbar(contour_plot_1,ax=ax_HA_NA,orientation='horizontal',ticks = [0,numpy.amax(histogram_list[0])],anchor=(0.5,1.4))
        cbar_2 = matplotlib.pyplot.colorbar(contour_plot_2,ax=ax_HA_M2,orientation='horizontal',ticks = [0,numpy.amax(histogram_list[1])],anchor=(0.5,1.4))
        ax_HA_NA.set_aspect('equal')
        ax_HA_M2.set_aspect('equal')
        matplotlib.pyplot.figtext(0.25,0.21,'frequency')
        matplotlib.pyplot.figtext(0.67,0.21,'frequency')
        fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/lipid_protein_proximity_2d_histo/' + outfile_name, dpi = 300)
    #call the plotting function for the different types of replicates (currently working with CHOL only):
    for replicate_string,histo_list,filename in zip(['sim35','sim36','sim37', 'sim38'],[sim_35_list_2d_histograms, sim_36_list_2d_histograms, sim_37_list_2d_histograms, sim_38_list_2d_histograms],['sim35_small_test.png','sim36_small_test.png','sim37_small_test.png','sim38_small_test.png']):
        replicate_contour_plotter('CHOL',replicate_string,histo_list,filename,'120')

def lipid_protein_closest_approach_analysis(queue_IPC,universe_object):
    '''This function should be similar to protein_proximity_analysis_v2() above, but instead of counting the number of proteins local to a given type of lipid, it should simply parse the closest COG-COG approach between a given lipid and a given protein type. The end objective is to produce a nice heat map for each replicate, and possibly to take heat map snapshots over time to see how things are progressing, but perhaps overall maps are fine to start with and directly compare the different replicates.'''
    #---------------------------------------------------------------------------------------
    #write a general function for producing a distance matrix between the residues of a given lipid type and a given protein type, and assessing the closest approach therein:
    def lipid_residue_protein_closest_approach(lipid_type, protein_type):
        '''For each of the residues of a given input lipid type and protein type, assess the COG-COG distance matrix and return a numpy array with the closest approach (A) between each lipid residue and the protein of the given type. 
        protein_type can be: 'HA', 'NA', or 'M2'
        lipid_type can be 'POPS', 'DOPX', 'DOPE', 'CHOL', 'PPCH', 'FORS'
        '''
        #the number of CG particles in the protein is dictated by protein_type, as is the number of proteins
        if protein_type == 'HA':
            protein_particles = 3621
            num_proteins = 80
            initial_topology_index = 1
            #perform a single MDA selection to pull out all of the HA coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 1:289680').coordinates()
        elif protein_type == 'NA':
            protein_particles = 4094
            num_proteins = 12
            initial_topology_index = 289681
            #perform a single MDA selection to pull out all of the NA coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 289681:338808').coordinates()
        elif protein_type == 'M2':
            protein_particles = 372
            num_proteins = 15
            initial_topology_index = 338809
            #perform a single MDA selection to pull out all of the M2 coordinates in the system:
            protein_coord_array = universe_object.selectAtoms('bynum 338809:344388').coordinates()

        #as was done with the proteins, I want to obtain a numpy array of coordinates for the all the lipids of the specified lipid type
        specified_lipid_selection_string = 'resname {lipid_residue_name}'.format(lipid_residue_name = lipid_type) 
        lipid_selection = universe_object.selectAtoms(specified_lipid_selection_string)
        lipid_coord_array = lipid_selection.coordinates()
        #I'm also going to want to know the number of particles than are within a single lipid residue for striding purposes:
        num_lipids = float(len(lipid_selection.resids()))
        particles_per_lipid = float(lipid_coord_array.shape[0]) / num_lipids  

        #now, I want a function that can 'stride' through the protein and lipid numpy coordinate arrays to average coordinate sets into the centers of geometry for individual proteins or lipid residues into new numpy arrays
        def produce_COG_arrays(particles_per_molecule, coordinate_array, num_molecules, molecule_type):
            '''Given a numpy coordinate_array of shape (num_particles, 3), use the particles_per_molecule to slice through the numpy array, average the Cartesian coordinates for each molecule/assembly to produce the COG of that assembly, and return an ordered numpy array of the COGs with shape (num_molecules, 3).
            molecule_type can be 'HA','NA','M2', or 'lipid'; the point is to only use TMDs for the protein COG calculations.'''
            COG_list = [] #append molecule/unit COG coords here, in topological order
            remaining_molecules = num_molecules 
            start_row_coefficient = 0
            def TMD_particle_selector(input_array):
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

            while remaining_molecules > 0:
                start_row = start_row_coefficient * particles_per_molecule
                end_row = ((start_row_coefficient + 1) * particles_per_molecule) 
                current_molecule_coordinate_array = coordinate_array[start_row:end_row,...]
                #deal with TMD selection if applicate (lipid cases will simply return the same array)
                current_molecule_coordinate_array = TMD_particle_selector(current_molecule_coordinate_array)
                #assert current_molecule_coordinate_array.shape == (particles_per_molecule,3), "current_molecule_coordinate_array has incorrect shape"
                #average together the above set of Cartesian coordinates for this molecule to obtain its COG:
                current_molecule_COG_array = numpy.average(current_molecule_coordinate_array,axis = 0)
                assert current_molecule_COG_array.shape == (3,), "current_molecule_COG_array has incorrect shape"
                #append the COG for the current molecule to the list and then increment/decrement varibles and continue loop:
                COG_list.append(current_molecule_COG_array)
                start_row_coefficient += 1 #pattern for row limits is 0:n , n:2n , 2n:3n , etc.
                remaining_molecules -= 1
            #convert the COG_list to a numpy array and verify its shape:
            COG_array = numpy.array(COG_list)
            assert COG_array.shape == (num_molecules, 3), "COG_array has incorrect shape."
            #finally return the numpy array containing the ordered set of per-molecule or per-assembly center of geometry values:
            return COG_array
    
        #use the above function to obtain two ordered per-molecule or per-assembly COG coordinate arrays:
        lipid_COG_coordinate_array = produce_COG_arrays(particles_per_lipid, lipid_coord_array, num_lipids, 'lipid')
        protein_COG_coordinate_array = produce_COG_arrays(protein_particles, protein_coord_array, num_proteins, protein_type)

        #I want to calculate the distance between each lipid of the specified type and each protein of the specified type
        import scipy, scipy.spatial.distance
        distance_matrix = scipy.spatial.distance.cdist(lipid_COG_coordinate_array,protein_COG_coordinate_array,'euclidean')
        assert distance_matrix.shape == (lipid_COG_coordinate_array.shape[0],protein_COG_coordinate_array.shape[0]), "distance_matrix has incorrect shape."
        #save distance_matrix to a text file for debugging purposes:
        #numpy.savetxt('distance_matrix_debug.txt',distance_matrix,newline ="\n")
        #each row of distance_matrix should correspond to the full set of COG-COG distances between the lipid residue (row) and each of the proteins of the given type (columns)

        closest_approach_array = distance_matrix.min(axis=1) #should return an ordered numpy array with the per-row (i.e., per-lipid) closest COG-COG approaches to the given protein type (A)
        assert closest_approach_array.shape == (lipid_COG_coordinate_array.shape[0],), "closest_approach_array has incorrect shape."
        #return the per-lipid closest approach array:
        return closest_approach_array
    #---------------------------------------------------------------------------------------

    range_list = [[0,300],[0,300]] #between 0 and 300 A range values for each dimension
    list_closest_approach_arrays= []
    for ts in universe_object.trajectory[::100]: #parse every Nth frame
        if ts.frame == 1: #for the first frame only
            for protein_type in ['HA','NA','M2']:
                list_closest_approach_arrays.append(lipid_residue_protein_closest_approach('CHOL',protein_type))
            #now, I want to generate a 2D histogram for HA relative to NA, and another 2D histogram for HA relative to M2:
            HA_NA_histogram, HA_edges_1, NA_edges = numpy.histogram2d(list_closest_approach_arrays[0],list_closest_approach_arrays[1],bins=301,range=range_list)
            HA_M2_histogram, HA_edges_2, M2_edges = numpy.histogram2d(list_closest_approach_arrays[0],list_closest_approach_arrays[2],bins=301,range=range_list)
            print multiprocessing.current_process().name, 'frame: ', ts.frame
        else:
            list_closest_approach_arrays = [] #reset for reuse 
            for protein_type in ['HA','NA','M2']:
                list_closest_approach_arrays.append(lipid_residue_protein_closest_approach('CHOL',protein_type))
            HA_NA_histogram_current, HA_edges_current, NA_edges_current = numpy.histogram2d(list_closest_approach_arrays[0],list_closest_approach_arrays[1],bins=301,range=range_list)
            HA_M2_histogram_current, HA_edges_2_current, M2_edges_current = numpy.histogram2d(list_closest_approach_arrays[0],list_closest_approach_arrays[2],bins=301,range=range_list)
            #I think it would not be good if the bins / edge list values changed during the trajectory (as summing histogram frequencies over the length of a trajectory would then become invalid), so check for this:
            assert numpy.subtract(HA_edges_current,HA_edges_1).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(NA_edges_current,NA_edges).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(HA_edges_2_current,HA_edges_2).sum() == 0, "histogram edge inconsistency"
            assert numpy.subtract(M2_edges_current,M2_edges).sum() == 0, "histogram edge inconsistency"
            #if the edge list assertions pass, then combine histograms:
            HA_NA_histogram = numpy.add(HA_NA_histogram,HA_NA_histogram_current)
            HA_M2_histogram = numpy.add(HA_M2_histogram,HA_M2_histogram_current)
            print multiprocessing.current_process().name, 'frame: ', ts.frame
        

    
    queue_IPC.put([[HA_NA_histogram,HA_edges_current, NA_edges_current],[HA_M2_histogram,HA_edges_2_current, M2_edges_current]]) #return the object to the IPC queue
    print 'object returned to the queue for:', multiprocessing.current_process().name #debugging


def capture_plot_lipid_protein_closest_approach_analysis(Q2,Q3,Q4,Q5,from_pickle = 'no'):
    '''Capture, modify, and plot the data from the matching multicore analysis function above. As the matching function only works sensibly in the presence of viral proteins, at present sims 35,36,37, and 38 apply, but sim 33 does not (so use queue objects 2,3,4,5).'''
    print 'start plotting code' #for debugging purposes
    #if not plotting from pickle files then generate objects directly from the multiprocessing queues and then pickle for safety/future use:
    if from_pickle == 'no':
        #capture the lists of 2d histograms from the IPC queue objects, with format: [[HA_NA_histogram,HA_edges_current, NA_edges_current],[HA_M2_histogram,HA_edges_2_current, M2_edges_current]]
        sim_35_list_2d_histograms= Q2.get()
        sim_36_list_2d_histograms= Q3.get()
        sim_37_list_2d_histograms= Q4.get()
        sim_38_list_2d_histograms= Q5.get()
        #pickle the lists to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim_35_list_2d_histograms,open('sim_35_list_2d_histograms_closest_contacts.p','wb'))
        pickle.dump(sim_36_list_2d_histograms,open('sim_36_list_2d_histograms_closest_contacts.p','wb'))
        pickle.dump(sim_37_list_2d_histograms,open('sim_37_list_2d_histograms_closest_contacts.p','wb'))
        pickle.dump(sim_38_list_2d_histograms,open('sim_38_list_2d_histograms_closest_contacts.p','wb'))
        print 'pickling complete' #debugging
    elif from_pickle == 'yes': #for this case I simply retrieve the lists from the pickled files
        list_of_lists= [pickle.load(open(x,'rb')) for x in ['sim_35_list_2d_histograms_closest_contacts.p','sim_36_list_2d_histograms_closest_contacts.p','sim_37_list_2d_histograms_closest_contacts.p','sim_38_list_2d_histograms_closest_contacts.p']]
        sim_35_list_2d_histograms, sim_36_list_2d_histograms, sim_37_list_2d_histograms, sim_38_list_2d_histograms = list_of_lists
    #with the histogram data now available, whether by direct retrieval from queue objects or by retrieval from pickled objects, go ahead and write some test plotting code to visualize the data:
    #write a general function that can be called to plot the data for each replicate:
    import matplotlib.pyplot
    def replicate_contour_plotter(lipid_name,replicate_name,histogram_object_list,outfile_name):
        '''replicate_name can be 'sim35','sim36','sim37', or 'sim38' 
        lipid_type can be 'POPS', 'DOPX', 'DOPE', 'CHOL', 'PPCH', 'FORS' '''
        fig = matplotlib.pyplot.figure()
        #testing: plot lipid species protein-proximity data for HA/NA and HA/M2 side-by-side
        ax_HA_NA = fig.add_subplot(121) #HA_NA data in first plot column 
        ax_HA_NA.set_xlabel('Closest COG approach to NA ($\AA$)')
        ax_HA_NA.set_ylabel('Closest COG approach to HA ($\AA$)')
        #matplotlib.pyplot.title(50 * ' ' + '{lipid_string} protein closest COG approach ({replicate_string} test case)'.format(replicate_string = replicate_name, lipid_string = lipid_name))
        ax_HA_M2 = fig.add_subplot(122) #HA_M2 data in second plot column 
        ax_HA_M2.set_xlabel('Closest COG approach to M2 ($\AA$)')
        #ax_HA_M2.set_ylabel('Closest COG approach to HA ($\AA$)')
        histogram_list = [histogram_object_list[0][0],histogram_object_list[1][0]]
        #obtain the extents from the queue/pickled objects as well:
        xedges_HA_NA = histogram_object_list[0][1]
        yedges_HA_NA = histogram_object_list[0][2]
        xedges_HA_M2 = histogram_object_list[1][1]
        yedges_HA_M2 = histogram_object_list[1][2]
        extent_list = [[yedges_HA_NA[0],yedges_HA_NA[-1],xedges_HA_NA[0],xedges_HA_NA[-1]],[yedges_HA_M2[0],yedges_HA_M2[-1],xedges_HA_M2[0],xedges_HA_M2[-1]]]
        contour_plot_1 = ax_HA_NA.contourf(histogram_list[0],extent=extent_list[0],origin='lower',levels = numpy.linspace(0,numpy.amax(histogram_list[0])))
        contour_plot_2 = ax_HA_M2.contourf(histogram_list[1],extent=extent_list[1],origin='lower',levels = numpy.linspace(0,numpy.amax(histogram_list[1])))
        cbar_1 = matplotlib.pyplot.colorbar(contour_plot_1,ax=ax_HA_NA,orientation='horizontal',ticks = [0,numpy.amax(histogram_list[0])],anchor=(0.5,1.3))
        cbar_2 = matplotlib.pyplot.colorbar(contour_plot_2,ax=ax_HA_M2,orientation='horizontal',ticks = [0,numpy.amax(histogram_list[1])],anchor=(0.5,1.3))
        ax_HA_NA.set_aspect('equal')
        ax_HA_M2.set_aspect('equal')
        matplotlib.pyplot.figtext(0.25,0.20,'frequency')
        matplotlib.pyplot.figtext(0.67,0.20,'frequency')
        fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/lipid_protein_closest_approach/' + outfile_name, dpi = 300)
    #call the plotting function for the different types of replicates (currently working with CHOL only):
    for replicate_string,histo_list,filename in zip(['sim35','sim36','sim37', 'sim38'],[sim_35_list_2d_histograms, sim_36_list_2d_histograms, sim_37_list_2d_histograms, sim_38_list_2d_histograms],['sim35_small_test.png','sim36_small_test.png','sim37_small_test.png','sim38_small_test.png']):
        replicate_contour_plotter('CHOL',replicate_string,histo_list,filename)

def euclidean_displacement_assessment(queue_IPC,universe_object,lipid_species):
    '''This function should assess the euclidean displacement of a given lipid species by tracing out the (x,y,z) movement of representative members that are: 1) the most mobile; 2) the closest to the average mobility; and 3) the least mobile. The main reason I'm writing this function is because I'm concerned that the lipids in my vesicle / virion simulations are not moving around very much, especially in the translational sense (see the 4D heat maps I've generated with concentric densities for each species). If I do find that there are at least some reaonably-mobile lipids, it may then be sensible to think about calculating the diffusion coefficient, perhaps in a fashion similar to equation 4 in Falck et al. (2004) Biophys J. 87: 1076-91, but save that for another function.'''
    import numpy.linalg
    def produce_lipid_species_selection(lipid_species_name):
        '''Simply produce the MDA selection object and calculate the num residues BEFORE trajectory iteration for a given lipid residue type, as this will get dynamically updated in the trajectory. It would be less efficient to re-select with each frame. lipid_species_name can be: DOPX, DOPE, POPS, CHOL, PPCH, FORS (input as a string). Returns a tuple (MDA_selection, num_residues)'''
        MDA_selection_string = 'resname {lipid_name}'.format(lipid_name = lipid_species_name)
        MDA_selection = universe_object.selectAtoms(MDA_selection_string)
        num_residues = MDA_selection.numberOfResidues()
        return (MDA_selection, num_residues)

    #define a series of utility functions that can be used during trajectory parsing:
    def produce_COG_array_current_frame(MDA_selection, num_residues):
        '''Return a numpy array of COG coordinates for the requested lipid species. The array should have shape (N_lipid_molecules,3).'''
        #print 'num_residues:', num_residues #debugging
        #produce a list of residue coordinate arrays by splitting the coordinate data according to the number of residues:
        list_lipid_residue_coord_arrays = numpy.split(MDA_selection.coordinates(),num_residues)
        #print "list_lipid_residue_coord_arrays.shape", numpy.array(list_lipid_residue_coord_arrays).shape #debugging
        #produce a numpy array containing the centroids of each lipid residue:
        array_lipid_residue_centroids = numpy.average(numpy.array(list_lipid_residue_coord_arrays),axis=1)
        #print "array_lipid_residue_centroids.shape", array_lipid_residue_centroids.shape #debugging
        assert array_lipid_residue_centroids.shape == (num_residues,3), "array_lipid_residue_centroids does not have the expected shape" 
        return array_lipid_residue_centroids 

    def produce_euclidean_displacement_array_current_frame(previous_array_lipid_residue_centroids,current_array_lipid_residue_centroids,num_residues):
        '''Return a numpy array containing the euclidean displacement of each lipid residue COG for lipid residues of the specified type compared to the previously assessed frame.'''
        #print "previous_array_lipid_residue_centroids.shape", previous_array_lipid_residue_centroids.shape #debugging
        #print "current_array_lipid_residue_centroids.shape", current_array_lipid_residue_centroids.shape #debugging
        array_lipid_residue_COG_euclidean_displacements = numpy.apply_along_axis(numpy.linalg.norm,1,previous_array_lipid_residue_centroids - current_array_lipid_residue_centroids)
        #print array_lipid_residue_COG_euclidean_displacements #debugging
        #print "array_lipid_residue_COG_euclidean_displacements.shape" , array_lipid_residue_COG_euclidean_displacements.shape 
        assert array_lipid_residue_COG_euclidean_displacements.shape == (num_residues,), "array_lipid_residue_COG_euclidean_displacements does not have the expected shape." #basically, it should have this shape because there should be a scalar displacement value for each residue centroid
        return array_lipid_residue_COG_euclidean_displacements

    def maintain_running_total_euclidean_displacement(current_euclidean_displacement_array,running_total_euclidean_displacement_array):
        '''Add in the most recent euclidean displacement for the current lipid species to the running total array.'''
        running_total_euclidean_displacement_array = numpy.add(running_total_euclidean_displacement_array,current_euclidean_displacement_array)
        return running_total_euclidean_displacement_array 

    
    def assess_final_euclidean_displacement(final_euclidean_displacement_array):
        '''This function will perform a variety of analyses on the final euclidean displacement array for a given lipid residue type and return the resulting objects in a tuple.'''
        #the max value in the final displacement array corresponds to the most mobile lipid residue of this type, so I'll want the total max displacement value and the Python index of this value (i.e, 0 is first residue of this type, etc.)
        max_euclidean_displacement_value = numpy.amax(final_euclidean_displacement_array)
        most_mobile_residue_index = numpy.where(final_euclidean_displacement_array == max_euclidean_displacement_value)
        #likewise for the minimum displacement / least mobile lipid residue of this type:
        min_euclidean_displacement_value = numpy.amin(final_euclidean_displacement_array)
        least_mobile_residue_index = numpy.where(final_euclidean_displacement_array == min_euclidean_displacement_value)
        #the average mobility will be a little different, as I'll want to find the residue with the closest (but probably not identical) mobility to the avg
        avg_euclidean_displacement_value = numpy.average(final_euclidean_displacement_array)
        #to find the representative avg mobility lipid residue, subtract the average from the final array and select the index of the minimum absolute value
        absolute_avg_difference_array = numpy.absolute(numpy.subtract(final_euclidean_displacement_array,avg_euclidean_displacement_value))
        avg_mobility_residue_index = numpy.where(absolute_avg_difference_array == numpy.amin(absolute_avg_difference_array))
        #I might as well histogram the euclidean mobility for this lipid type as well:
        hist, bin_edges = numpy.histogram(final_euclidean_displacement_array,bins=20)
        #now return all the values in a tuple:
        return (max_euclidean_displacement_value, most_mobile_residue_index, min_euclidean_displacement_value, least_mobile_residue_index, avg_euclidean_displacement_value, avg_mobility_residue_index, hist, bin_edges)


    def store_relevant_Cartesian_coords(array_lipid_residue_centroids,most_mobile_residue_index,least_mobile_residue_index,avg_mobility_residue_index,cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile):
        '''The overall code will iterate once through a given trajectory to assess which residues have the most mobility, least mobility, and closest-to-average mobility, and then a second iteration through the same trajectory with this function will be used to capture the per-frame Cartesian coordinates. Note that I'll only want to capture coordinates for the 3 lipid residue COGs of interest, so that I can plot the displacement trace for the representative cases.'''
        #index the relevant set of Cartesian coordinates (for the representative lipid residues) and concatenate their values into the respective arrays:
        cartesian_storage_array_most_mobile = numpy.concatenate((cartesian_storage_array_most_mobile,array_lipid_residue_centroids[most_mobile_residue_index,...][0,...]))
        cartesian_storage_array_least_mobile = numpy.concatenate((cartesian_storage_array_least_mobile,array_lipid_residue_centroids[least_mobile_residue_index,...][0,...]))
        cartesian_storage_array_avg_mobile = numpy.concatenate((cartesian_storage_array_avg_mobile,array_lipid_residue_centroids[avg_mobility_residue_index,...][0,...]))
        #verify that there are only 3 columns of data as the traj iteration progresses:
        #print "cartesian_storage_array_most_mobile.shape", cartesian_storage_array_most_mobile.shape #debugging
        assert cartesian_storage_array_most_mobile.shape[1] == 3, "cartesian_storage_array_most_mobile has incorrect shape"
        assert cartesian_storage_array_least_mobile.shape[1] == 3, "cartesian_storage_array_least_mobile has incorrect shape"
        assert cartesian_storage_array_avg_mobile.shape[1] == 3, "cartesian_storage_array_avg_mobile has incorrect shape"
        #return the concatenated arrays in a tuple:
        return (cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile)


    #-----------------------
    #first pass through trajectory, for assessing the mobility of a population of lipid residues of the specified type:
    #can produce the MDA selection and count the residues of the specified lipid type BEFORE trajectory iteration starts:
    MDA_lipid_selection, num_lipid_residues = produce_lipid_species_selection(lipid_species)
    running_total_euclid_array = numpy.zeros(num_lipid_residues) #initialize to keep track of running total euclidean displacements below
    for ts in universe_object.trajectory: #every N frames for now
        #start off by producing an array of lipid residue centroid coordinates for the current frame:
        #for the first frame, I have to do something a bit different, as there's no previous set of centroids:
        if ts.frame == 1:
            prev_array_lipid_res_centroids = produce_COG_array_current_frame(MDA_lipid_selection, num_lipid_residues)
            print multiprocessing.current_process().name, 'first pass frame:', ts.frame #monitor progress on terminal
        elif ts.frame%100 == 0: #skip 100 frames this way as it's faster in MDA at the moment than using [::] notation ;for all subsequent frames there will be a 'previous' reference set of centroids for displacement calculation
            current_array_lipid_res_centroids = produce_COG_array_current_frame(MDA_lipid_selection, num_lipid_residues)
            #compare current and previous lipid centroid arrays to obtain euclidean displacements
            array_lipid_residue_COG_euclidean_displacements = produce_euclidean_displacement_array_current_frame(prev_array_lipid_res_centroids,current_array_lipid_res_centroids,num_lipid_residues)
            #maintain a running total of the euclidean displacements for each lipid residue (for the second parsed frame, this will simply add to an array of zeros)
            running_total_euclid_array = maintain_running_total_euclidean_displacement(array_lipid_residue_COG_euclidean_displacements,running_total_euclid_array)
            #now, before continuing the loop, I'll want to set the current array of centroids to previous status, as they will be carried over in the next frame
            prev_array_lipid_res_centroids = current_array_lipid_res_centroids
            print multiprocessing.current_process().name, 'first pass frame:', ts.frame #monitor progress on terminal
    #after trajectory iteration completes for this replicate, analyse the running total euclidean distance array:
    parameter_tuple = assess_final_euclidean_displacement(running_total_euclid_array)
    #where parameter_tuple has format: max_euclidean_displacement_value, most_mobile_residue_index, min_euclidean_displacement_value, least_mobile_residue_index, avg_euclidean_displacement_value, avg_mobility_residue_index, hist, bin_edges 

    #-----------------------
    #second pass through trajectory, to capture the Cartesian coordinates for the representative lipid residues:
    for ts in universe_object.trajectory: #every N frames for now
        if ts.frame == 1: #first frame only
            #start off by producing the overall array of lipid residue centroids:
            array_lipid_res_centroids = produce_COG_array_current_frame(MDA_lipid_selection, num_lipid_residues)
            #then initialize the relevant Cartesian coordinate arrays with the values from the first frame:
            cartesian_storage_array_most_mobile = array_lipid_res_centroids[parameter_tuple[1],...][0,...]
            cartesian_storage_array_least_mobile = array_lipid_res_centroids[parameter_tuple[3],...][0,...]
            cartesian_storage_array_avg_mobile = array_lipid_res_centroids[parameter_tuple[5],...][0,...]
            print multiprocessing.current_process().name, 'second pass frame:', ts.frame #monitor progress on terminal
        elif ts.frame%100 == 0: #for all subsequent frames (every 100th frame)
            #start off by producing the overall array of lipid residue centroids:
            array_lipid_res_centroids = produce_COG_array_current_frame(MDA_lipid_selection, num_lipid_residues)
            #concatenate the relevant COG coords to the appropriate array:
            cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile = store_relevant_Cartesian_coords(array_lipid_res_centroids,parameter_tuple[1],parameter_tuple[3],parameter_tuple[5],cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile)        
            print multiprocessing.current_process().name, 'second pass frame:', ts.frame #monitor progress on terminal


    #-----------------------
    #now that the two trajectory loopings have completed, go ahead and combine objects as needed and return them to the queue so they can be pickled for safety / future use:
    queue_IPC.put([parameter_tuple,cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile])


def capture_plot_euclidean_displacement_assessment(Q1,Q2,Q3,Q4,Q5,from_pickle = 'no',lipid_name = 'CHOL'):
    '''Retrieve data from the queue objects produced by the similarly-named analysis function and pickle / plot it.'''
    if from_pickle == 'no':
        #format for data lists: parameter_tuple,cartesian_storage_array_most_mobile,cartesian_storage_array_least_mobile,cartesian_storage_array_avg_mobile
        #where parameter_tuple has format: max_euclidean_displacement_value, most_mobile_residue_index, min_euclidean_displacement_value, least_mobile_residue_index, avg_euclidean_displacement_value, avg_mobility_residue_index, hist, bin_edges 
        sim33_data_list = Q1.get()
        sim35_data_list = Q2.get()
        sim36_data_list = Q3.get()
        sim37_data_list = Q4.get()
        sim38_data_list = Q5.get()
        #pickle the lists to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim33_data_list,open('sim33_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'wb'))
        pickle.dump(sim35_data_list,open('sim35_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'wb'))
        pickle.dump(sim36_data_list,open('sim36_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'wb'))
        pickle.dump(sim37_data_list,open('sim37_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'wb'))
        pickle.dump(sim38_data_list,open('sim38_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'wb'))
        print 'pickling complete' #debugging
    else: #load the data from the stored pickle files
        sim33_data_list, sim35_data_list,sim36_data_list,sim37_data_list, sim38_data_list = [pickle.load(open(x,'rb')) for x in ['sim33_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'sim35_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'sim36_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'sim37_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name),'sim38_euclidean_displacement_data_{lipid_name}.p'.format(lipid_name=lipid_name)]]
    #with the data available either directly from the parsed trajectories or from saved pickle files, I can now focus on the plotting code
    import matplotlib.pyplot,matplotlib.pylab
    matplotlib.pylab.rcParams['ytick.major.pad']='20'
    def euclidean_panel_plotter(data_list,sim_name,lipid_name):
        '''A general function for plotting the panel analysis layouts from a given replicate and for a given lipid type.'''
        #start by plotting the 3 representative Cartesian-tracking plots--------------------
        fig = matplotlib.pyplot.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax_least_mobile_lipid_residue_trace = fig.add_subplot(231,projection='3d',aspect="equal") 
        ax_avg_mobile_lipid_residue_trace = fig.add_subplot(232,projection='3d',aspect="equal") 
        ax_most_mobile_lipid_residue_trace = fig.add_subplot(233,projection='3d',aspect="equal") 
        figure_least_title='Least mobile {lipid_resname} residue {simulation_name}\n total displacement: {displacement}'.format(displacement=str(round(data_list[0][2],1)),lipid_resname=lipid_name,simulation_name=sim_name)
        figure_avg_title='Avg mobility {lipid_resname} residue {simulation_name}\n total displacement: {displacement}'.format(displacement=str(round(data_list[0][-4],1)),lipid_resname=lipid_name,simulation_name=sim_name)
        figure_most_title='Most mobile {lipid_resname} residue {simulation_name}\n total displacement: {displacement}'.format(displacement=str(round(data_list[0][0],1)),lipid_resname=lipid_name,simulation_name=sim_name)
        ax_least_mobile_lipid_residue_trace.scatter(data_list[2][:,0].ravel(),data_list[2][:,1].ravel(),data_list[2][:,2].ravel(),zdir='z',c='r',marker='d')
        ax_avg_mobile_lipid_residue_trace.scatter(data_list[3][:,0].ravel(),data_list[3][:,1].ravel(),data_list[3][:,2].ravel(),zdir='z',c='r',marker='d')
        ax_most_mobile_lipid_residue_trace.scatter(data_list[1][:,0].ravel(),data_list[1][:,1].ravel(),data_list[1][:,2].ravel(),zdir='z',c='r',marker='d')
        for axis in [ax_most_mobile_lipid_residue_trace,ax_avg_mobile_lipid_residue_trace,ax_least_mobile_lipid_residue_trace]:
             axis.set_xlim3d(400,1100)
             axis.set_ylim3d(400,1100)
             axis.set_zlim3d(-50,650)
             axis.set_xticks(numpy.arange(400,1200,200))
             axis.set_yticks(numpy.arange(500,1300,200)) #trying to shift these so they are placed properly; matplotlib currently does a bad job with 3D plot space management
             #axis.set_yticks(numpy.arange(400,1200,200)) #trying to shift these so they are placed properly; matplotlib currently does a bad job with 3D plot space management
             axis.set_yticklabels(numpy.arange(400,1200,200))
             axis.set_zticks(numpy.arange(0,700,200))
             axis.set_xlabel('\n\nx' + r'($\AA$)')
             axis.set_ylabel('\n\ny' + r'($\AA$)')
             axis.set_zlabel('\n\nz' + r'($\AA$)')

        #adjust title placements manually
        ax_least_mobile_lipid_residue_trace.text(50,750,880,figure_least_title + r' $\AA$')
        ax_avg_mobile_lipid_residue_trace.text(150,750,900,figure_avg_title + r' $\AA$')
        ax_most_mobile_lipid_residue_trace.text(320,750,940,figure_most_title + r' $\AA$')
        #end Cartesian tracking plot code--------------------------------

        #now start the histogram plot code---------------------------------
        ax_histo = fig.add_subplot(212)
        center = (data_list[0][-1][:-1]+data_list[0][-1][1:])/2
        width = 0.95*(data_list[0][-1][1]-data_list[0][-1][0])
        ax_histo.bar(center,data_list[0][-2],align='center',width=width)
        ax_histo.set_xlabel('Cumulative Euclidean Displacement ' + r'($\AA$)')
        ax_histo.set_ylabel('# {resname} residues'.format(resname=lipid_name))
        ax_histo.set_xlim(1500,6000)
        ax_histo.set_ylim(0,2000)
        #ax_histo.annotate('least', xy=(data_list[0][2], 0),  xycoords='data', xytext=(data_list[0][2] - width/2.0, -numpy.max(data_list[0][-2]/4.0)), textcoords='data', arrowprops=dict(arrowstyle="->"))
        ax_histo.annotate('least', xy=(data_list[0][2], 0),  xycoords='data', xytext=(data_list[0][2] - width/2.0, -numpy.max(data_list[0][-2]/4.0) -290.0), textcoords='data', arrowprops=dict(arrowstyle="->"))
        #ax_histo.annotate('avg', xy=(data_list[0][-4],0),  xycoords='data', xytext=(data_list[0][-4] - width/2.0, -numpy.max(data_list[0][-2]/4.0)),textcoords='data', arrowprops=dict(arrowstyle="->"))
        ax_histo.annotate('avg', xy=(data_list[0][-4],0),  xycoords='data', xytext=(data_list[0][-4] - width/2.0, -numpy.max(data_list[0][-2]/4.0) -290.0),textcoords='data', arrowprops=dict(arrowstyle="->"))
        #ax_histo.annotate('most', xy=(data_list[0][0],0),  xycoords='data', xytext=(data_list[0][0] - width/2.0, -numpy.max(data_list[0][-2]/4.0)), textcoords='data', arrowprops=dict(arrowstyle="->"))
        ax_histo.annotate('most', xy=(data_list[0][0],0),  xycoords='data', xytext=(data_list[0][0] - width/2.0, -numpy.max(data_list[0][-2]/4.0) -290.0), textcoords='data', arrowprops=dict(arrowstyle="->"))
        #end histogram plot code------------------------------------------

        #fig.set_size_inches(18,12) #careful with sizing to maintain roughly equal plot edges; 3D aspect ratio not properly built-in for matplotlib
        fig.set_size_inches(11,7.33) #careful with sizing to maintain roughly equal plot edges; 3D aspect ratio not properly built-in for matplotlib
        fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/euclidean_displacement/{sim_name}_euclidean_displacement_{lipid_type}.png'.format(sim_name=sim_name,lipid_type=lipid_name), dpi = 300)

    #now call the code for the set of replicates and lipid type you want to plot:
    for data,sim in zip([sim33_data_list, sim35_data_list,sim36_data_list,sim37_data_list, sim38_data_list],['sim33','sim35','sim36','sim37','sim38']):
        euclidean_panel_plotter(data,sim,lipid_name)
    print 'plotting complete for', lipid_name


def sphericity_tracking(queue_IPC,universe_object):
    '''Tracking the sphericity of the equilibrating vesicle or the production virion simulations using lipid molecules only (protein spikes would obviously completely alter the sphericity measure). Based on the sphericity description available at http://en.wikipedia.org/wiki/Sphericity , which cites the original primary source as Wadell (1935). See additional notes / details in the validated prototype module I import and use below.'''
    import sphericity_tracking_prototype #the module I used to design / test / validate most of the utility functions for sphericity analysis
    import scipy,scipy.spatial
    #before trajectory iteration, make the MDA lipid selection (shouldn't matter if some replicates are missing a few lipids--they simply won't get selected):
    #there seem to be two lipid residues 'flating around' in the vesicle sims; these were deleted for the virion sims, but I want to selectively ignore these residues when parsing vesicle sphericity:
    num_protein_residues = universe_object.selectAtoms('protein').numberOfResidues()
    if num_protein_residues > 0: #virion simulation, so just use the normal selection
        MDA_lipid_selection = universe_object.selectAtoms('resname PPCH or resname FORS or resname POPS or resname DOPE or resname DOPX or resname CHOL')
    else: #vesicle simulation, so ignore the two offending floater residues, which could throw off sphericity otherwise
        MDA_lipid_selection = universe_object.selectAtoms('(resname PPCH or resname FORS or resname POPS or resname DOPE or resname DOPX or resname CHOL) and not (resid 8896 or resid 18204)')

    #now define a list where sphericity values can be stored in order and start iterating through the trajectory:
    list_sphericity_values_in_trajectory = []
    for ts in universe_object.trajectory:
        if ts.frame == 1 or ts.frame%10 == 0: #only use first or every 10th frame
            #get the full set of lipid coordinates in the current frame:
            lipid_coordinates_current_frame = MDA_lipid_selection.coordinates()
            #produce a scipy convex hull object for the full set of lipid coordinates in this frame:
            convex_hull_scipy_object = scipy.spatial.ConvexHull(lipid_coordinates_current_frame)
            #use my prototyped function to calculate the sphericity using the convex hull object:
            sphericity_value_this_frame = sphericity_tracking_prototype.sphericity_proc(convex_hull_scipy_object)
            #append the sphericity value to the list for this trajectory:
            list_sphericity_values_in_trajectory.append(sphericity_value_this_frame)
            #track progress at the terminal:
            print multiprocessing.current_process().name, 'frame:', ts.frame , '[sphericity:', sphericity_value_this_frame, ']'
        else: continue
    #after iterating through the trajectory on this particular core, return the array of sphericity values to the appropriate queue object, where they will then be pickled for analysis / plotting later:
    queue_IPC.put(numpy.array(list_sphericity_values_in_trajectory))


def capture_plot_sphericity_tracking(Q1,Q2,Q3,Q4,Q5,Q6,from_pickle = 'no'):
    '''Capture / pickle /plot the values returned to queue objects from the similarly named parent analysis function.'''
    if from_pickle == 'no':
        sim33_sphericity_tracking_array = Q1.get()
        sim35_sphericity_tracking_array = Q2.get()
        sim36_sphericity_tracking_array = Q3.get()
        sim37_sphericity_tracking_array = Q4.get()
        sim38_sphericity_tracking_array = Q5.get()
        sim39_sphericity_tracking_array = Q6.get()
        #pickle the arrays to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
        pickle.dump(sim33_sphericity_tracking_array,open('sim33_sphericity.p','wb'))
        pickle.dump(sim35_sphericity_tracking_array,open('sim35_sphericity.p','wb'))
        pickle.dump(sim36_sphericity_tracking_array,open('sim36_sphericity.p','wb'))
        pickle.dump(sim37_sphericity_tracking_array,open('sim37_sphericity.p','wb'))
        pickle.dump(sim38_sphericity_tracking_array,open('sim38_sphericity.p','wb'))
        pickle.dump(sim39_sphericity_tracking_array,open('sim39_sphericity.p','wb'))
        print 'pickling complete' #debugging
    else: #load the data from the stored pickle files
        sim33_sphericity_tracking_array,sim35_sphericity_tracking_array,sim36_sphericity_tracking_array,sim37_sphericity_tracking_array,sim38_sphericity_tracking_array,sim39_sphericity_tracking_array = [pickle.load(open(x,'rb')) for x in ['sim33_sphericity.p','sim35_sphericity.p','sim36_sphericity.p','sim37_sphericity.p','sim38_sphericity.p','sim39_sphericity.p']]
    #plotting code here when ready:
    import matplotlib, matplotlib.pyplot
    matplotlib.rc('xtick', labelsize=8)
    matplotlib.rc('ytick', labelsize=8)
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    colour_dict = {'vesicle':'#000000','virion\n (295 K)':'#006600','virion restrained\n (295 K)':'#FF0000','virion\n(323 K)':'#00FF00','FORS\nvirion (295 K)':'#000066','FORS\nvirion (323 K)':'#00CCFF'}
    for sphericity_tracking_array,sim_name in zip([sim33_sphericity_tracking_array,sim35_sphericity_tracking_array,sim36_sphericity_tracking_array,sim37_sphericity_tracking_array,sim38_sphericity_tracking_array,sim39_sphericity_tracking_array],['vesicle','virion\n (295 K)','virion restrained\n (295 K)','virion\n(323 K)','FORS\nvirion (295 K)','FORS\nvirion (323 K)']):
        colour_value = colour_dict[sim_name]
        if sim_name == 'vesicle': #no need to adjust the frame numbers / time for the equilibration and vesicle sim
            ax.plot(numpy.arange(sphericity_tracking_array.size)/1000.0,sphericity_tracking_array,label=sim_name,linewidth=2,c=colour_value)
        else: #for all subsequent, virion-based, simulations I want to add 300 ns (0.3 us) to make it clear that I'm starting from the equilibrated vesicle end point
            ax.plot(numpy.arange(sphericity_tracking_array.size)/1000.0 + 0.3,sphericity_tracking_array,label=sim_name,linewidth=2,c=colour_value)
    ax.set_xlabel('Time ($\mu$s)',fontsize=14)
    ax.set_ylabel('$\Psi$',fontsize=18)
    ax.set_ylim((0.982,0.998)) 
    #ax.set_xlim((-0.02,0.5)) #a first plot emphasing the early / equilibration stages of the trajectories
    #ax.legend(bbox_to_anchor=(-.14, 1.04, 1.25, .102), loc=3,
               #ncol=6, mode="expand", borderaxespad=0.)
    matplotlib.rc('legend',**{'fontsize':8,'columnspacing':0})
    ax.legend(loc=4,ncol=2)
    matplotlib.pyplot.subplots_adjust(left=0.25, right=0.97, top=0.90, bottom=0.15)
    #fig.set_size_inches(11,4.0) #for smaller plots to be used in SysL meeting slide
    fig.set_size_inches(4.0,3.5) #for smaller plots to be used in SysL meeting slide
    #fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/sphericity_tracking/sphericity_tracking_equilibration.png',dpi=300)
    ax.set_xlim((-0.2,5.3)) #a second plot over the full course of the trajectories
    fig.savefig('/sansom/sc2/bioc1009/Dropbox/Tyler/flu_MS_figure_work/Fig_2_parts/Figure_2B_300_DPI.png',dpi=300)
    #for seminar:
    #matplotlib.pyplot.subplots_adjust(left=0.15, right=0.97, top=0.90, bottom=0.15)
    #fig.set_size_inches(11,4.0)
    #fig.savefig('/sansom/sc2/bioc1009/python_scripts/matplotlib_scripts/output_plots/sphericity_tracking/sphericity_seminar_full.png',dpi=300)
        

def mean_square_displacement(queue_IPC,universe_object):
    '''Calculate the MSD for the various trajectories as part of the lipid diffusion analysis.'''
    import scipy, scipy.spatial.distance
    #MSD_window_size_frames = 500 #I think 5000 is about 1/10th of a 5 us trajectory (i.e., about every 500 ns)
    #incorporate a new check for the window size (frame skip value) list because Danny's trajectory data skips every 10th frame
    if universe_object.selectAtoms('resname DUPC').numberOfResidues() != 0:
        MSD_window_size_list = [1,3,5,10,25,50,100,200,300,400,500] #Danny's sim has SKIP 10 filtering so use 1/10th of normal frame values; dt = 1 ns to 500 ns
    else:
        MSD_window_size_list = [10,30,50,100,250,500,1000,2000,3000,4000,5000] #My sims; dt = 1 ns to 500 ns

    num_protein_residues = universe_object.selectAtoms('protein').numberOfResidues()
    if num_protein_residues > 0: #virion simulation, so just use the normal selection
        selection_string_all_lipids = 'resname DOPE or resname DOPX or resname PPCH or resname CHOL or resname POPS or resname FORS or resname DPPC or resname DUPC'
    else: #vesicle simulation, so ignore the two offending floater residues, which could throw off diffusion otherwise
        selection_string_all_lipids = '(resname PPCH or resname FORS or resname POPS or resname DOPE or resname DOPX or resname CHOL) and not (resid 8896 or resid 18204)'

    MDA_selection_all_lipids = universe_object.selectAtoms(selection_string_all_lipids) #updates automatically during trajectory looping
    #write a function to produce a numpy array of the centroids of each lipid residue
    def centroid_array_production(current_MDA_selection):
        residue_list = current_MDA_selection.residues
        centroid_list = [residue.centroid() for residue in residue_list]
        return numpy.array(centroid_list)
    #write a function to pull out the coordinates of the appropriate headgroup particle of each lipid residue
    def headgroup_coordinate_array_production(current_MDA_selection):
        headgroup_particle_selection = current_MDA_selection.selectAtoms('((resname POPS or resname DOPX or resname DOPE or resname PPCH or resname DPPC or resname DUPC) and name PO4) or (resname CHOL and name ROH) or (resname FORS and name AM2)')
        headgroup_particle_coordinate_array = headgroup_particle_selection.coordinates()
        return headgroup_particle_coordinate_array

    dict_MSD_values_this_replicate = {'MSD_value_list':[],'frame_skip_value_list':[],'MSD_value_list_headgroups':[]}
    for MSD_window_size_frames in MSD_window_size_list:
        list_per_window_average_displacements = []
        list_per_window_average_displacements_headgroups = []
        counter = 0
        for ts in universe_object.trajectory[::MSD_window_size_frames]:
            if counter == 0: #first parsed frame
                previous_frame_lipid_residue_centroid_array = centroid_array_production(MDA_selection_all_lipids)
                previous_frame_lipid_residue_headgroup_coord_array = headgroup_coordinate_array_production(MDA_selection_all_lipids)
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            else: #all subsequent frames
                current_frame_lipid_residue_centroid_array = centroid_array_production(MDA_selection_all_lipids)
                current_frame_lipid_residue_headgroup_array = headgroup_coordinate_array_production(MDA_selection_all_lipids)
                assert current_frame_lipid_residue_centroid_array.shape == current_frame_lipid_residue_headgroup_array.shape, "Lipid residue centroid & heagroup arrays should have the same shape."
                #calculate the MSD here:
                #allocate a numpy array for storage of [r(t) - r(0)]**2 displacements
                #square_displacement_storage_array = numpy.zeros(current_frame_lipid_residue_centroid_array.shape[0])
                #a bit awkward, but just iterate through both position arrays and calculate the pairwise distances
                #improve square displacement algorithm by using array operations as much as possible:
                delta_array = previous_frame_lipid_residue_centroid_array - current_frame_lipid_residue_centroid_array
                delta_array_headgroups = previous_frame_lipid_residue_headgroup_coord_array - current_frame_lipid_residue_headgroup_array
                square_delta_array = numpy.square(delta_array)
                square_delta_array_headgroups = numpy.square(delta_array_headgroups)
                sum_squares_delta_array = numpy.sum(square_delta_array,axis=1)
                sum_squares_delta_array_headgroups = numpy.sum(square_delta_array_headgroups,axis=1)
                MSD_all_lipids_this_frame = numpy.average(sum_squares_delta_array)
                MSD_all_lipids_this_frame_headgroups = numpy.average(sum_squares_delta_array_headgroups)
                list_per_window_average_displacements.append(MSD_all_lipids_this_frame)
                list_per_window_average_displacements_headgroups.append(MSD_all_lipids_this_frame_headgroups)
                #reset the value of the 'previous' array as you go along:
                previous_frame_lipid_residue_centroid_array = current_frame_lipid_residue_centroid_array
                previous_frame_lipid_residue_headgroup_coord_array = current_frame_lipid_residue_headgroup_array
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            counter += 1
        average_displacement_value_current_dt = numpy.average(numpy.array(list_per_window_average_displacements))
        average_displacement_value_current_dt_headgroup = numpy.average(numpy.array(list_per_window_average_displacements_headgroups))
        dict_MSD_values_this_replicate['MSD_value_list'].append(average_displacement_value_current_dt)
        dict_MSD_values_this_replicate['MSD_value_list_headgroups'].append(average_displacement_value_current_dt_headgroup)
        dict_MSD_values_this_replicate['frame_skip_value_list'].append(MSD_window_size_frames)
    queue_IPC.put(dict_MSD_values_this_replicate)

def capture_plot_MSD(Q1,Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no',sim_danny_only='no'):
    '''For now I'm just pickling the mean square displacement data, but will eventually want to plot it.'''
    if sim_danny_only == 'yes': #only analyzing Danny's simulation data so just pickle it & exit (will deal with plotting from full set of pickle files later)
        sim_Danny_MSD_tracking_dictionary = Q7.get()
        pickle.dump(sim_Danny_MSD_tracking_dictionary,open('sim_Danny_MSD.p','wb'))
    else: #regular control flow if not exclusively dealing with Danny's virion sim
        if from_pickle == 'no':
            sim33_MSD_tracking_dictionary = Q1.get()
            sim35_MSD_tracking_dictionary = Q2.get()
            sim36_MSD_tracking_dictionary = Q3.get()
            sim37_MSD_tracking_dictionary = Q4.get()
            sim38_MSD_tracking_dictionary = Q5.get()
            sim39_MSD_tracking_dictionary = Q6.get()
            #pickle the dictionaries to be safe and so that I can easily adjust plots without re-parsing all the trajectories:
            pickle.dump(sim33_MSD_tracking_dictionary,open('sim33_MSD.p','wb'))
            pickle.dump(sim35_MSD_tracking_dictionary,open('sim35_MSD.p','wb'))
            pickle.dump(sim36_MSD_tracking_dictionary,open('sim36_MSD.p','wb'))
            pickle.dump(sim37_MSD_tracking_dictionary,open('sim37_MSD.p','wb'))
            pickle.dump(sim38_MSD_tracking_dictionary,open('sim38_MSD.p','wb'))
            pickle.dump(sim39_MSD_tracking_dictionary,open('sim39_MSD.p','wb'))
            print 'pickling complete' #debugging
        else: #load the data from the stored pickle files
            sim33_MSD_tracking_dictionary,sim35_MSD_tracking_dictionary,sim36_MSD_tracking_dictionary,sim37_MSD_tracking_dictionary,sim38_MSD_tracking_dictionary,sim39_MSD_tracking_dictionary,sim_Danny_MSD_tracking_dictionary = [pickle.load(open(x,'rb')) for x in ['sim33_MSD.p','sim35_MSD.p','sim36_MSD.p','sim37_MSD.p','sim38_MSD.p','sim39_MSD.p','sim_Danny_MSD.p']]
        #now plot / analyze the data
        import matplotlib, matplotlib.pyplot
        fig = matplotlib.pyplot.figure()
        subplot_number = 241
        list_D_values = []
        for MSD_tracking_dictionary,sim_name in zip([sim33_MSD_tracking_dictionary,sim35_MSD_tracking_dictionary,sim36_MSD_tracking_dictionary,sim37_MSD_tracking_dictionary,sim38_MSD_tracking_dictionary,sim39_MSD_tracking_dictionary,sim_Danny_MSD_tracking_dictionary],['sim33 (323 K / 295 K vesicle equil)','sim35 (295 K virion)','sim36 (295 K virion -- proteins restrained)','sim37 (323 K virion)','sim38 (FORS 295 K)','sim39 (FORS 323 K)','sim Danny (40% CHOL)']):
            ax = fig.add_subplot(subplot_number)
            ax.set_title(sim_name)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('MSD ($\AA^2$)')
            #ax.set_xticks([1,3,5,10,25,50,100,200,300,400,500])
            ax.set_xticks([100,200,300,400,500])
            ax.set_yticks(numpy.arange(0,6000,1000))
            ax.set_xlim((0,500))
            ax.set_ylim((0,6000))
            array_MSD_values = numpy.array(MSD_tracking_dictionary['MSD_value_list'])
            array_MSD_values_headgroups = numpy.array(MSD_tracking_dictionary['MSD_value_list_headgroups'])
            array_frame_values = numpy.array(MSD_tracking_dictionary['frame_skip_value_list'])
            if not sim_name == 'sim Danny (40% CHOL)':
                array_ns_time_points = array_frame_values / 10.0
            else: #for Danny's virion sim, the values are already equivalent to ns (because of the SKIP 10 processing)
                array_ns_time_points = array_frame_values 
            slope, intercept = numpy.polyfit(array_ns_time_points,array_MSD_values,1) #assuming that these are actually lines and not higher order polynomial trends
            headgroup_slope, headgroup_intercept = numpy.polyfit(array_ns_time_points,array_MSD_values_headgroups,1) #assuming that these are actually lines and not higher order polynomial trends
            coefficients = slope, intercept
            headgroup_coefficients = headgroup_slope, headgroup_intercept
            polynomial = numpy.poly1d(coefficients)
            polynomial_headgroup = numpy.poly1d(headgroup_coefficients)
            x = numpy.arange(0,600,100)
            ys = polynomial(x)
            ys_headgroup = polynomial_headgroup(x)
            ax.plot(x,ys,color='black') #best fit line
            ax.plot(x,ys_headgroup,color='red') #best fit line
            #divide the slope (A ^ 2 / ns ) by 10 ^ 7 to convert to standard cm ^ 2 / s units:
            diffusion_constant = slope / (10 ** 7)
            diffusion_constant_headgroup = headgroup_slope / (10 ** 7)
            list_D_values.extend([diffusion_constant,diffusion_constant_headgroup])
            diffusion_constant = '%.3e' % diffusion_constant #convert to scientific notation
            diffusion_constant_headgroup = '%.3e' % diffusion_constant_headgroup #convert to scientific notation
            diffusion_constant_label_string = 'D$_{{centroids}}$ = {diff_constant} cm$^2$ / s'.format(diff_constant = diffusion_constant)
            diffusion_constant_label_string_headgroup = 'D$_{{headgroups}}$ = {diff_constant} cm$^2$ / s'.format(diff_constant = diffusion_constant_headgroup)
            text_1_y_position = 5500
            text_2_y_position = 4900
            if sim_name == 'sim Danny (40% CHOL)': #line break here to position text on this particular plot
                diffusion_constant_label_string_headgroup = 'D$_{{headgroups}}$ = \n{diff_constant} cm$^2$ / s'.format(diff_constant = diffusion_constant_headgroup)
                text_2_y_position = 4600
            ax.text(30, text_1_y_position, diffusion_constant_label_string, fontsize=12,bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'),color='black')
            ax.text(30, text_2_y_position, diffusion_constant_label_string_headgroup, fontsize=12,bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'),color='red')
            ax.scatter(array_ns_time_points,array_MSD_values,color='black')
            ax.scatter(array_ns_time_points,array_MSD_values_headgroups,color='red')
            subplot_number += 1
        #now I want an 8th panel plot that compares our diffusion constant values to those previously reported for biological systems in the literature
        ax = fig.add_subplot(subplot_number)
        ax.set_title('Literature Comparison')
        ax.set_xlabel('% CHOL')
        ax.set_xlim((-1,70))
        ax.set_ylabel('D (cm$^2$ / s)')
        CHOL_values_our_sims = numpy.array([52] * 12 + [40,40]) #52% for all of my virion sims & 40% for Danny's
        #print CHOL_values_our_sims
        D_values_our_sims = numpy.array(list_D_values)
        #print D_values_our_sims
        ax.scatter(CHOL_values_our_sims,D_values_our_sims,label='our results',color='black')
        ax.scatter(numpy.array([52,52]),numpy.array([7.0,35.0])/ (10**8),label='flu experiment',color='red') #solid-state NMR results between 45 C & 17 C from Polozov et al. Nature Chemical Biology 4: 248-55. [assume 52% CHOL because this is influenza proper (or at least its lipids)]
        ax.scatter(numpy.array([0,25]),numpy.array([21.0,1.0]) / (10**8), label='experiment',color='green') #NMR results for a quaternary bilayer between 303K -- 333K at a variety of CHOL % values (only using those values corresponding to D extremities here) from Fillipov et al. (2004) Biophys J. 86: 891-6.
        ax.scatter(numpy.array([10.0,10.0,20.0,20.0,33.0,33.0]),numpy.array([4.9,0.105,5.15,0.255,5.1,0.795])/(10**8),color='green') #FCS Lo & Ld phase results from Kahya et al. (2003) JBC 278: 28109-28115.
        ax.scatter(numpy.array([0,15,30,45,60]),numpy.array([3.0,1.8,1.1,0.6,0.25])/(10**8),color='green') #FCS on GUVs by Korlach et al. (1999) PNAS 96: 8461-6.
        ax.scatter(numpy.array([0]),numpy.array([8.5])/(10**7),label='published CG',color='blue') #CG-MD by Goose & Sansom (2013) PLoS Comput Biol 9(4): e1003033. 
        ax.scatter(numpy.array([0,4.7,12.5,20.3,29.7,50.0]),numpy.array([10.9,7.9,3.9,1.5,0.9,0.2])/(10**8),label='published AT',color='orange') #AT-MD by Falck et al. (2004) Biophys J. 87: 1076-91.
        ax.scatter(numpy.array([0]),numpy.array([1.2])/(10**7),color='orange') #Lindahl and Edholm (2001) J Chem Phys 115: 4938.
        ax.scatter(numpy.array([0]),numpy.array([3.0])/(10**7),color='orange') #Essmann and Berkowitz (1999) Biophys J 76: 2081-9.
        ax.scatter(numpy.array([0,0,5,5,15,15,20,20,25,25,27.5,27.5,35,35,37.5,37.5]),numpy.array([5.1,27.0,3.0,24.0,1.5,20.0,1.1,17.5,0.9,16.0,0.9,15.1,0.9,12.0,0.9,12.0])/(10**8),color='green') #values at T extremes from Flippov et al. (2003) Biophys J. 84: 3079-86.
        ax.scatter(numpy.array([0,40,0,40,0,40,0,40,0,40]),numpy.array([3.43,2.29,3.55,2.38,5.53,2.46,5.68,2.61,6.08,2.75])/(10**8),color='orange') #MD by Rabinovich et al. (2007) Biochemistry (Moscow) 1: 343-357.
        ax.scatter(numpy.array([0,40,0,40,0,40,0,40]),numpy.array([7.99,4.40,10.3,5.78,12.8,6.30,14.1,8.55])/(10**8),color='green') #NMR by Rabinovich et al. (2007) Biochemistry (Moscow) 1: 343-357.
        ax.scatter(numpy.array([0,0,0,0]),numpy.array([1.67,1.27,1.30,1.46])/(10**7),color='orange') #AT-MD by Flenner et al. (2009) Physical Review E 79: 
        ax.scatter(numpy.array([0,0,0,0]),numpy.array([1.25,1.15,1.07,1.06])/(10**7),color='brown',label='theory') #theory by Flenner et al. (2009) Physical Review E 79: 
        ax.scatter(numpy.array([0,0]),numpy.array([5.8,5.7])/(10**7),color='blue') #CG-MD by Stachura and Kneller (2013) Molecular Simulation 40: 
        ax.scatter(numpy.array([0]),numpy.array([1.8])/(10**7),color='orange') #AT-MD by Stachura and Kneller (2013) Molecular Simulation 40: 
        ax.scatter(numpy.zeros(4),numpy.array([5.62,4.85,6.15,4.53])/(10**7),color='blue') #CG-MD by Karjiban et al. (2013) (Hindawi) Journal of Chemistry Article ID 931051, 6 pages, 2013. doi:10.1155/2013/931051
        ax.scatter(numpy.zeros(4),numpy.array([1.0,4.0,0.005,0.04])/10**7,color='blue') #CG-MD by Marrink et al. (2005) Chemistry and Physics of Lipids 135: 223-244.
        ax.scatter(numpy.array([15,15]),numpy.array([3.5,0.5])/(10**8),color='orange') #AT MD by Hakobyan and Heuer (2013) J Phys Chem B 117: 3841-51.
        ax.scatter(numpy.array([15]),numpy.array([60.0])/(10**8),color='blue') #CG MD by Hakobyan and Heuer (2013) J Phys Chem B 117: 3841-51.
        ax.scatter(numpy.ones(6) * 20,numpy.array([6.8,4.5,9.0,10.5,5.0,5.2])/(10**8),color='blue') #CG-MD by Li and Gorfe (2013) PLoS ONE 8(7): e71018. doi:10.1371/journal.pone.0071018
        ax.scatter(numpy.array([0,0,10,10,10,10,25,25,25,25,25,25,25,25]),numpy.array([4.15,4.06,0.33,3.83,2.45,4.61,0.16,2.16,1.02,2.13,0.16,4.10,0.36,1.93])/(10**7),color='blue') #CG-MD by Apajalahti et al. (2010) Faraday Discussions 144: 411-430.
        ax.scatter(numpy.array([25]),numpy.array([0.08])/(10**7),color='orange') #AT-MD by Apajalahti et al. (2010) Faraday Discussions 144: 411-430.
        
        ax.legend(loc='upper right',scatterpoints=1,ncol=1,prop={'size':9})
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        fig.set_size_inches(20,9) 
        matplotlib.pyplot.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)
        fig.savefig('MSD_vs_time_influenza.png', dpi = 300)


def mean_square_displacement_by_species(queue_IPC,universe_object):
    '''Calculate the MSD for the various trajectories as part of the lipid diffusion analysis. This is a modified version of the function designed to break the analysis down according to lipid species.'''
    import scipy, scipy.spatial.distance
    #MSD_window_size_frames = 500 #I think 5000 is about 1/10th of a 5 us trajectory (i.e., about every 500 ns)
    #incorporate a new check for the window size (frame skip value) list because Danny's trajectory data skips every 10th frame
    if universe_object.selectAtoms('resname DUPC').numberOfResidues() != 0:
        MSD_window_size_list = [1,3,5,10,25,50,100,200,300,400,500] #Danny's sim has SKIP 10 filtering so use 1/10th of normal frame values; dt = 1 ns to 500 ns
    else:
        MSD_window_size_list = [10,30,50,100,250,500,1000,2000,3000,4000,5000] #My sims; dt = 1 ns to 500 ns

    num_protein_residues = universe_object.selectAtoms('protein').numberOfResidues()
    if num_protein_residues > 0: #virion simulation, so just use the normal selection
        #selection_string_all_lipids = 'resname DOPE or resname DOPX or resname PPCH or resname CHOL or resname POPS or resname FORS or resname DPPC or resname DUPC'
        #list_lipid_selection_strings = ['resname DOPE','resname DOPX', 'resname PPCH','resname CHOL','resname POPS', 'resname FORS','resname DPPC','resname DUPC']
        lipid_selection_dictionary = {'DOPE':{'residue_selection':'resname DOPE','residue_headgroup_selection':'resname DOPE and name PO4'},'DOPX':{'residue_selection':'resname DOPX','residue_headgroup_selection':'resname DOPX and name PO4'},'PPCH':{'residue_selection':'resname PPCH','residue_headgroup_selection':'resname PPCH and name PO4'},'CHOL':{'residue_selection':'resname CHOL','residue_headgroup_selection':'resname CHOL and name ROH'},'POPS':{'residue_selection':'resname POPS','residue_headgroup_selection':'resname POPS and name PO4'},'FORS':{'residue_selection':'resname FORS','residue_headgroup_selection':'resname FORS and name AM2'},'DPPC':{'residue_selection':'resname DPPC','residue_headgroup_selection':'resname DPPC and name PO4'},'DUPC':{'residue_selection':'resname DUPC','residue_headgroup_selection':'resname DUPC and name PO4'}}
    else: #vesicle simulation, so ignore the two offending floater residues, which could throw off diffusion otherwise
        #selection_string_all_lipids = '(resname PPCH or resname FORS or resname POPS or resname DOPE or resname DOPX or resname CHOL) and not (resid 8896 or resid 18204)'
        #list_lipid_selection_strings = ['resname DOPE and not (resid 8896 or resid 18204)','resname DOPX and not (resid 8896 or resid 18204)', 'resname PPCH and not (resid 8896 or resid 18204)','resname CHOL and not (resid 8896 or resid 18204)','resname POPS and not (resid 8896 or resid 18204)', 'resname FORS and not (resid 8896 or resid 18204)','resname DPPC and not (resid 8896 or resid 18204)','resname DUPC and not (resid 8896 or resid 18204)']
        lipid_selection_dictionary = {'DOPE':{'residue_selection':'resname DOPE and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname DOPE and not (resid 8896 or resid 18204) and name PO4'},'DOPX':{'residue_selection':'resname DOPX and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname DOPX and not (resid 8896 or resid 18204) and name PO4'},'PPCH':{'residue_selection':'resname PPCH and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname PPCH and not (resid 8896 or resid 18204) and name PO4'},'CHOL':{'residue_selection':'resname CHOL and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname CHOL and not (resid 8896 or resid 18204) and name ROH'},'POPS':{'residue_selection':'resname POPS and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname POPS and not (resid 8896 or resid 18204) and name PO4'},'FORS':{'residue_selection':'resname FORS and not (resid 8896 or resid 18204)','residue_headgroup_selection':'resname FORS and not (resid 8896 or resid 18204) and name AM2'},'DPPC':{'residue_selection':'resname DPPC','residue_headgroup_selection':'resname DPPC and name PO4'},'DUPC':{'residue_selection':'resname DUPC','residue_headgroup_selection':'resname DUPC and name PO4'}}
        #lipid_selection_dictionary = {'DOPE':'resname DOPE and not (resid 8896 or resid 18204)','DOPX':'resname DOPX and not (resid 8896 or resid 18204)','PPCH':'resname PPCH and not (resid 8896 or resid 18204)','CHOL':'resname CHOL and not (resid 8896 or resid 18204)','POPS':'resname POPS and not (resid 8896 or resid 18204)','FORS':'resname FORS and not (resid 8896 or resid 18204)','DPPC':'resname DPPC','DUPC':'resname DUPC'}

    #MDA_selection_all_lipids = universe_object.selectAtoms(selection_string_all_lipids) #updates automatically during trajectory looping
    #MDA_residue_selection_list = [universe_object.selectAtoms(residue_selection) for residue_selection in list_lipid_selection_strings] #updates automatically during trajectory looping
    MDA_residue_selection_dictionary = {}
    for lipid_name,selection_subdictionary in lipid_selection_dictionary.iteritems():
        current_MDA_residue_selection = universe_object.selectAtoms(selection_subdictionary['residue_selection'])
        current_MDA_residue_headgroup_selection = universe_object.selectAtoms(selection_subdictionary['residue_headgroup_selection'])
        #print 'lipid name', lipid_name
        #print 'MDA sel num residues:', current_MDA_residue_selection.numberOfResidues()
        if current_MDA_residue_selection.numberOfResidues() > 0:
            MDA_residue_selection_dictionary[lipid_name] = {'residue_selection':current_MDA_residue_selection,'residue_headgroup_selection':current_MDA_residue_headgroup_selection}
    print 'MDA_residue_selection_dictionary.keys():',MDA_residue_selection_dictionary.keys()
    #so MDA_residue_selection_dictionary should have key = lipid name, value = subdictionary with 'residue_selection' & 'residue_headgroup_selection' entries

    #write a function to produce a numpy array of the centroids of each lipid residue (sublisted by residue type now)
    def centroid_array_production(current_MDA_selection_dictionary):
        dictionary_centroid_arrays = {}
        for lipid_name,selection_subdictionary in current_MDA_selection_dictionary.iteritems():
            #print lipid_name
            current_list_residues = selection_subdictionary['residue_selection'].residues
            #print 'num residues:', len(current_list_residues)
            current_list_residue_centroids = [residue.centroid() for residue in current_list_residues]
            dictionary_centroid_arrays[lipid_name] = numpy.array(current_list_residue_centroids)
        return dictionary_centroid_arrays
    #write a function to pull out the coordinates of the appropriate headgroup particle of each lipid residue (sublisted by residue type now)
    def headgroup_coordinate_array_production(current_MDA_selection_dictionary):
        headgroup_coordinate_array_dictionary = {}
        for lipid_name,selection_subdictionary in current_MDA_selection_dictionary.iteritems():
            #print lipid_name
            current_list_residues = selection_subdictionary['residue_selection'].residues
            current_headgroup_selection = selection_subdictionary['residue_headgroup_selection']
            #print 'current_headgroup_selection', current_headgroup_selection
            current_headgroup_coordinate_array = current_headgroup_selection.coordinates()
            headgroup_coordinate_array_dictionary[lipid_name] = current_headgroup_coordinate_array
        return headgroup_coordinate_array_dictionary

    dict_MSD_values_this_replicate = {'MSD_value_dict':{},'MSD_std_dict':{},'frame_skip_value_list':[],'MSD_value_dict_headgroups':{},'MSD_std_dict_headgroups':{}} #for overall storage of MSD average / standard deviation values for this replicate
    for MSD_window_size_frames in MSD_window_size_list:
        list_per_window_average_displacements = []
        list_per_window_average_displacements_headgroups = []
        counter = 0
        trajectory_striding_dictionary = {} #store values in a dict as you stride through trajectory
        for ts in universe_object.trajectory[::MSD_window_size_frames]:
            if counter == 0: #first parsed frame
                previous_frame_lipid_residue_centroid_array_dictionary = centroid_array_production(MDA_residue_selection_dictionary)
                previous_frame_lipid_residue_headgroup_coord_array_dictionary = headgroup_coordinate_array_production(MDA_residue_selection_dictionary)
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            else: #all subsequent frames
                current_frame_lipid_residue_centroid_array_dictionary = centroid_array_production(MDA_residue_selection_dictionary)
                current_frame_lipid_residue_headgroup_coord_array_dictionary = headgroup_coordinate_array_production(MDA_residue_selection_dictionary)
                #assert current_frame_lipid_residue_centroid_array.shape == current_frame_lipid_residue_headgroup_array.shape, "Lipid residue centroid & heagroup arrays should have the same shape."
                #calculate the MSD here:
                #allocate a numpy array for storage of [r(t) - r(0)]**2 displacements
                #square_displacement_storage_array = numpy.zeros(current_frame_lipid_residue_centroid_array.shape[0])
                #a bit awkward, but just iterate through both position arrays and calculate the pairwise distances
                #improve square displacement algorithm by using array operations as much as possible:
                #dictionary adjustment performed up to here-----------------------------
                
                #start drafting the per-frame dict-based algorithm:
                for lipid_name in current_frame_lipid_residue_centroid_array_dictionary.keys():
                    if not lipid_name in trajectory_striding_dictionary.keys(): #create the appropriate entry if this lipid species hasn't been parsed yet for this replicate
                        trajectory_striding_dictionary[lipid_name] = {'MSD_value_list_centroids':[],'MSD_value_list_headgroups':[]}
                        current_delta_array_centroids = previous_frame_lipid_residue_centroid_array_dictionary[lipid_name] - current_frame_lipid_residue_centroid_array_dictionary[lipid_name]
                        current_delta_array_headgroups = previous_frame_lipid_residue_headgroup_coord_array_dictionary[lipid_name] - current_frame_lipid_residue_headgroup_coord_array_dictionary[lipid_name]
                        square_delta_array_centroids = numpy.square(current_delta_array_centroids)
                        square_delta_array_headgroups = numpy.square(current_delta_array_headgroups)
                        sum_squares_delta_array_centroids = numpy.sum(square_delta_array_centroids,axis=1)
                        sum_squares_delta_array_headgroups = numpy.sum(square_delta_array_headgroups,axis=1)
                        trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'].append(numpy.average(sum_squares_delta_array_centroids))
                        trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'].append(numpy.average(sum_squares_delta_array_headgroups))
                    else:
                        current_delta_array_centroids = previous_frame_lipid_residue_centroid_array_dictionary[lipid_name] - current_frame_lipid_residue_centroid_array_dictionary[lipid_name]
                        current_delta_array_headgroups = previous_frame_lipid_residue_headgroup_coord_array_dictionary[lipid_name] - current_frame_lipid_residue_headgroup_coord_array_dictionary[lipid_name]
                        square_delta_array_centroids = numpy.square(current_delta_array_centroids)
                        square_delta_array_headgroups = numpy.square(current_delta_array_headgroups)
                        sum_squares_delta_array_centroids = numpy.sum(square_delta_array_centroids,axis=1)
                        sum_squares_delta_array_headgroups = numpy.sum(square_delta_array_headgroups,axis=1)
                        trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'].append(numpy.average(sum_squares_delta_array_centroids))
                        trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'].append(numpy.average(sum_squares_delta_array_headgroups))

                #reset the value of the 'previous' array as you go along:
                previous_frame_lipid_residue_centroid_array_dictionary = current_frame_lipid_residue_centroid_array_dictionary
                previous_frame_lipid_residue_headgroup_coord_array_dictionary = current_frame_lipid_residue_headgroup_coord_array_dictionary
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            counter += 1
        for lipid_name,MSD_data_subdictionary in trajectory_striding_dictionary.iteritems():
            if not lipid_name in dict_MSD_values_this_replicate['MSD_value_dict'].keys(): #initialize lipid species subdictionaries as needed
                dict_MSD_values_this_replicate['MSD_value_dict'][lipid_name] = []
                dict_MSD_values_this_replicate['MSD_std_dict'][lipid_name] = []
                dict_MSD_values_this_replicate['MSD_value_dict_headgroups'][lipid_name] = []
                dict_MSD_values_this_replicate['MSD_std_dict_headgroups'][lipid_name] = []
                dict_MSD_values_this_replicate['MSD_value_dict'][lipid_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_std_dict'][lipid_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_value_dict_headgroups'][lipid_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'])))
                dict_MSD_values_this_replicate['MSD_std_dict_headgroups'][lipid_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'])))
            else:
                #print 'MSD value list centroids:',trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids']
                #print 'value to append:', numpy.average(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids']))
                #print 'list to append to:', dict_MSD_values_this_replicate['MSD_value_dict'][lipid_name]
                dict_MSD_values_this_replicate['MSD_value_dict'][lipid_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_std_dict'][lipid_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_value_dict_headgroups'][lipid_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'])))
                dict_MSD_values_this_replicate['MSD_std_dict_headgroups'][lipid_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[lipid_name]['MSD_value_list_headgroups'])))
        dict_MSD_values_this_replicate['frame_skip_value_list'].append(MSD_window_size_frames)
    queue_IPC.put(dict_MSD_values_this_replicate)

def capture_plot_MSD_by_lipid(Q1,Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no'):
    '''Capture (and eventually plot) the per-lipid MSD data.'''
    if from_pickle == 'no': #capture the data directly from the queue objects:
        sim33_MSD_tracking_dictionary_by_species = Q1.get()
        sim35_MSD_tracking_dictionary_by_species = Q2.get()
        sim36_MSD_tracking_dictionary_by_species = Q3.get()
        sim37_MSD_tracking_dictionary_by_species = Q4.get()
        sim38_MSD_tracking_dictionary_by_species = Q5.get()
        sim39_MSD_tracking_dictionary_by_species = Q6.get()
        sim_Danny_MSD_tracking_dictionary_by_species = Q7.get()
        #pickle for safety / future plotting
        pickle.dump(sim33_MSD_tracking_dictionary_by_species,open('sim33_MSD_by_species.p','wb'))
        pickle.dump(sim35_MSD_tracking_dictionary_by_species,open('sim35_MSD_by_species.p','wb'))
        pickle.dump(sim36_MSD_tracking_dictionary_by_species,open('sim36_MSD_by_species.p','wb'))
        pickle.dump(sim37_MSD_tracking_dictionary_by_species,open('sim37_MSD_by_species.p','wb'))
        pickle.dump(sim38_MSD_tracking_dictionary_by_species,open('sim38_MSD_by_species.p','wb'))
        pickle.dump(sim39_MSD_tracking_dictionary_by_species,open('sim39_MSD_by_species.p','wb'))
        pickle.dump(sim_Danny_MSD_tracking_dictionary_by_species,open('sim_Danny_MSD_by_species.p','wb'))
    else: #load the data from pickled storage
        sim33_MSD_tracking_dictionary_by_species,sim35_MSD_tracking_dictionary_by_species,sim36_MSD_tracking_dictionary_by_species,sim37_MSD_tracking_dictionary_by_species,sim38_MSD_tracking_dictionary_by_species,sim39_MSD_tracking_dictionary_by_species,sim_Danny_MSD_tracking_dictionary_by_species = [pickle.load(open(x,'rb')) for x in ['sim33_MSD_by_species.p','sim35_MSD_by_species.p','sim36_MSD_by_species.p','sim37_MSD_by_species.p','sim38_MSD_by_species.p','sim39_MSD_by_species.p','sim_Danny_MSD_by_species.p']]
        #(eventually) plot the data by lipid species
        import matplotlib.pyplot
        #adjust the value of index_window_filter as needed for dealing with anomalous diffusion****
        def register_diffusion_constant_plot_best_fit_line(output_dictionary, x_data_array, y_data_array, lipid_label,axis,index_window_filter):
            '''index_window_filter is for excluding values at short MSD time windows which may cause anomalous diffusion.'''
            x_data_array = x_data_array[index_window_filter:]
            y_data_array = y_data_array[index_window_filter:]
            slope, intercept = numpy.polyfit(x_data_array,y_data_array,1) #assuming that these are actually lines and not higher order polynomial trends
            #divide the slope (A ^ 2 / ns ) by 10 ^ 7 to convert to standard cm ^ 2 / s units:
            diffusion_constant = slope / (10 ** 7)
            #I also want to estimate the error in the diffusion constant (slope) values, and I'm going to mimic the approach of g_msd, which uses the difference between the slopes obtained for the two halves for the data:
            first_half_x_data, second_half_x_data = numpy.array_split(x_data_array,2)
            first_half_y_data, second_half_y_data = numpy.array_split(y_data_array,2)
            slope_first_half, intercept_first_half = numpy.polyfit(first_half_x_data,first_half_y_data,1) 
            slope_second_half, intercept_second_half = numpy.polyfit(second_half_x_data,second_half_y_data,1) 
            diffusion_constant_error_estimate = abs(slope_first_half - slope_second_half) / (10 ** 7)
            #put in the appropriate dictionary:
            output_dictionary[lipid_label_string] = {'diffusion_constant':diffusion_constant,'error_estimate':diffusion_constant_error_estimate}
            #plot the best linear fit to the data for this particular lipid species:
            coefficients = slope, intercept
            polynomial = numpy.poly1d(coefficients)
            x = numpy.arange(0,600,100)
            ys = polynomial(x)
            axis.plot(x,ys,color='black') #best fit line
        def diffusion_constant_bar_chart(diffusion_constant_dictionary,title,axis):
            x_positions = numpy.arange(0,len(diffusion_constant_dictionary))
            color_list = ['blue','red','purple','black','green','cyan']
            bars_plotted = 0
            for lipid_label,diffusion_dictionary in diffusion_constant_dictionary.iteritems():
                #diffusion_constant = '%.3e' % diffusion_constant #convert to scientific notation
                axis.bar(x_positions[bars_plotted],diffusion_dictionary['diffusion_constant'],label=lipid_label,color=color_list[bars_plotted],width=0.7,yerr=diffusion_dictionary['error_estimate'],ecolor='black')
                print lipid_label, diffusion_dictionary['diffusion_constant']
                axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axis.legend(bbox_to_anchor=(0.5,-0.05),loc='upper center',numpoints=1,prop={'size':6},ncol=6,fancybox = True)
                bars_plotted += 1
            axis.set_title(title)
            axis.set_ylabel('D (cm$^2$ / s)')
            axis.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
            axis.set_xlim(0,6)
            #axis.set_ylim(0,2e-6)
        fig = matplotlib.pyplot.figure()
        subplot_number = 1
        for MSD_tracking_dictionary,plot_name in zip([sim33_MSD_tracking_dictionary_by_species,sim35_MSD_tracking_dictionary_by_species,sim36_MSD_tracking_dictionary_by_species,sim37_MSD_tracking_dictionary_by_species,sim38_MSD_tracking_dictionary_by_species,sim39_MSD_tracking_dictionary_by_species,sim_Danny_MSD_tracking_dictionary_by_species],['vesicle 295/323 K', 'virion 295 K', 'protein-restrained virion 295 K','virion 323 K','FORS virion 295 K','FORS virion 323 K','Danny viron (40% CHOL) 310 K']):
            print plot_name
            headgroup_diffusion_constant_dictionary = {}
            centroid_diffusion_constant_dictionary = {}
            #print 'MSD tracking dictionary keys:', MSD_tracking_dictionary.keys()
            ax = fig.add_subplot(7,4,subplot_number)
            #print 'MSD value subdictionary keys:', MSD_tracking_dictionary['MSD_value_dict'].keys()
            #print 'MSD value subdictionary values:', MSD_tracking_dictionary['MSD_value_dict'].values()
            frame_window_array = numpy.array(MSD_tracking_dictionary['frame_skip_value_list'])
            if frame_window_array.max() == 500: #Danny's sim already skips 10 frames so convert directly to ns
                time_window_array_nanoseconds = frame_window_array
            else: #need to divide by 10 for other replicates
                time_window_array_nanoseconds = frame_window_array / 10.0
            for lipid_label_string, MSD_centroid_value_array in MSD_tracking_dictionary['MSD_value_dict'].iteritems():
                #print 'time window array:',time_window_array_nanoseconds
                #print 'MSD value array:',MSD_centroid_value_array
                #print 'MSD std dev array:',MSD_tracking_dictionary['MSD_std_dict'][lipid_label_string]
                ax.errorbar(time_window_array_nanoseconds,MSD_centroid_value_array,yerr=MSD_tracking_dictionary['MSD_std_dict'][lipid_label_string],label='{lipid_label} centroids'.format(lipid_label=lipid_label_string),ls='none',marker='o')
                register_diffusion_constant_plot_best_fit_line(output_dictionary=centroid_diffusion_constant_dictionary, x_data_array=time_window_array_nanoseconds, y_data_array=MSD_centroid_value_array, lipid_label =lipid_label_string,axis=ax,index_window_filter=0)
                #and a similar plot for headgroup-based MSD vs. time:
                ax2 = fig.add_subplot(7,4,subplot_number + 1)
                ax2.errorbar(time_window_array_nanoseconds,MSD_tracking_dictionary['MSD_value_dict_headgroups'][lipid_label_string],yerr=MSD_tracking_dictionary['MSD_std_dict_headgroups'][lipid_label_string],label='{lipid_label} headgroups'.format(lipid_label=lipid_label_string),ls='none',marker='o')
                #register_diffusion_constant_plot_best_fit_line(output_dictionary=headgroup_diffusion_constant_dictionary, x_data_array=time_window_array_nanoseconds, y_data_array=MSD_tracking_dictionary['MSD_value_dict_headgroups'][lipid_label_string], lipid_label =lipid_label_string,axis=ax2,index_window_filter=0)
            for ax in [ax, ax2]:
                ax.set_title(plot_name)
                ax.set_xlabel('Time (ns)')
                ax.set_ylabel('MSD ($\AA^2$)')
                ax.set_xlim(0,550)
                #ax.set_ylim(-2000,9000)
                ax.legend(loc=2,numpoints=1,prop={'size':6})

            #add two more columns for the D constant comparison bar charts
            ax3 = fig.add_subplot(7,4,subplot_number + 2)
            diffusion_constant_bar_chart(diffusion_constant_dictionary=centroid_diffusion_constant_dictionary,title='centroid D comparison',axis=ax3)
            ax4 = fig.add_subplot(7,4,subplot_number + 3)
            diffusion_constant_bar_chart(diffusion_constant_dictionary=headgroup_diffusion_constant_dictionary,title='headgroup D comparison',axis=ax4)
            subplot_number += 4
        fig.set_size_inches(32,40)
        #fig.savefig('MSD_analysis_by_lipid_species_flu.png',dpi=300) #use all data points
        #fig.savefig('MSD_analysis_by_lipid_species_flu_windows_over_50ns.png',dpi=300) #only use time windows > 50 ns
            
def plot_log_MSD_by_lipid_scaling_exponent():
    '''Use pickled per-lipid MSD data to plot log(MSD) vs. log(time) and determine the per-lipid scaling exponents (alpha values).'''
    sim33_MSD_tracking_dictionary_by_species,sim35_MSD_tracking_dictionary_by_species,sim36_MSD_tracking_dictionary_by_species,sim37_MSD_tracking_dictionary_by_species,sim38_MSD_tracking_dictionary_by_species,sim39_MSD_tracking_dictionary_by_species,sim_Danny_MSD_tracking_dictionary_by_species = [pickle.load(open(x,'rb')) for x in ['sim33_MSD_by_species.p','sim35_MSD_by_species.p','sim36_MSD_by_species.p','sim37_MSD_by_species.p','sim38_MSD_by_species.p','sim39_MSD_by_species.p','sim_Danny_MSD_by_species.p']]

    def calculate_scaling_exponent_plot_best_fit_line(log_x_data_array,log_y_data_array,axis,lipid_label):
        '''For calculation of scaling exponent (alpha) and plotting of the best fit line on the log plots.'''
        slope, intercept = numpy.polyfit(log_x_data_array,log_y_data_array,1) #assuming that these are actually lines and not higher order polynomial trends
        scaling_exponent = slope
        coefficients = slope, intercept
        polynomial = numpy.poly1d(coefficients)
        x = numpy.arange(0,3,0.1)
        ys = polynomial(x)
        axis.plot(x,ys,color='black') #best fit line
        #axis.text(0.5,0.5,'$\\alpha = {scaling_exponent}$ {lipid_label}'.format(scaling_exponent = scaling_exponent, lipid_label = protein_label))
        scaling_exponent_string = '%.2f' % scaling_exponent
        alpha_string = '$\\alpha = {scaling_exponent}$ {lipid_label}'.format(scaling_exponent = scaling_exponent_string, lipid_label = lipid_label)
        return alpha_string

    import matplotlib.pyplot
    colordict = {'DPPC':'black','DUPC':'pink','CHOL':'green','PPCH':'orange','DOPE':'purple','DOPX':'blue','POPS':'red','FORS':'brown'}
    fig = matplotlib.pyplot.figure()
    subplot_number = 1
    for MSD_tracking_dictionary,plot_name in zip([sim33_MSD_tracking_dictionary_by_species,sim35_MSD_tracking_dictionary_by_species,sim36_MSD_tracking_dictionary_by_species,sim37_MSD_tracking_dictionary_by_species,sim38_MSD_tracking_dictionary_by_species,sim39_MSD_tracking_dictionary_by_species,sim_Danny_MSD_tracking_dictionary_by_species],['vesicle 295/323 K', 'virion 295 K', 'protein-restrained virion 295 K','virion 323 K','FORS virion 295 K','FORS virion 323 K','virion (40% CHOL + DPPC / DUPC) 310 K']):
        print 'plot_name:',plot_name
        ax = fig.add_subplot(7,2,subplot_number)
        frame_window_array = numpy.array(MSD_tracking_dictionary['frame_skip_value_list'])
        if frame_window_array.max() == 500: #Danny's sim already skips 10 frames so convert directly to ns
            time_window_array_nanoseconds = frame_window_array
        else: #need to divide by 10 for other replicates
            time_window_array_nanoseconds = frame_window_array / 10.0
        for lipid_label_string, MSD_centroid_value_array in MSD_tracking_dictionary['MSD_value_dict'].iteritems():
            alpha_string = calculate_scaling_exponent_plot_best_fit_line(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_centroid_value_array),ax,lipid_label_string)
            ax.scatter(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_centroid_value_array),label='{alpha_string}'.format(alpha_string=alpha_string),marker='o',color=colordict[lipid_label_string])
            #and a similar plot for headgroup-based log(MSD) vs. log(time):
            ax2 = fig.add_subplot(7,2,subplot_number + 1)
            alpha_string = calculate_scaling_exponent_plot_best_fit_line(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_tracking_dictionary['MSD_value_dict_headgroups'][lipid_label_string]),ax2,lipid_label_string)
            ax2.scatter(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_tracking_dictionary['MSD_value_dict_headgroups'][lipid_label_string]),label='headgroup {alpha_string}'.format(alpha_string=alpha_string),marker='o',color=colordict[lipid_label_string])
        for ax in [ax, ax2]:
            ax.set_title(plot_name)
            ax.set_xlabel('log(t)')
            ax.set_ylabel('log(MSD)')
            ax.set_xlim(-0.1,3.0)
            ax.set_ylim(0,4.5)
            ax.legend(loc=2,scatterpoints=1,prop={'size':10})
        subplot_number += 2
    fig.set_size_inches(11,34)
    fig.savefig('log_MSD_scaling_exponents_per_lipid.png',dpi=300) 


    
def mean_square_displacement_proteins(queue_IPC,universe_object):
    '''MSD analysis code for analyzing the diffusion of the three different categories of protein in the flu virion constructs. Obviously, it only makes sense to use this code for the virion simulations & not the vesicle simulations without proteins.'''
    import scipy, scipy.spatial.distance
    #incorporate a check for the window size (frame skip value) list because Danny's trajectory data skips every 10th frame
    if universe_object.selectAtoms('resname DUPC').numberOfResidues() != 0:
        MSD_window_size_list = [1,3,5,10,25,50,100,200,300,400,500] #Danny's sim has SKIP 10 filtering so use 1/10th of normal frame values; dt = 1 ns to 500 ns
    else:
        MSD_window_size_list = [10,30,50,100,250,500,1000,2000,3000,4000,5000] #My sims; dt = 1 ns to 500 ns

    protein_selection_dictionary = {'all_HA_atoms_selection_string':'bynum 1:289680','all_NA_atoms_selection_string':'bynum 289681:338808','all_M2_atoms_selection_string':'bynum 338809:344388'}

    MDA_full_protein_group_selection_dictionary = {} #for storing protein groups by type, not individual proteins yet
    for selection_descriptor, selection_string in protein_selection_dictionary.iteritems():
        current_MDA_all_atom_selection_this_protein_type = universe_object.selectAtoms(selection_string)
        MDA_full_protein_group_selection_dictionary[selection_descriptor.replace('_string','')] = current_MDA_all_atom_selection_this_protein_type

    def centroid_array_production(current_MDA_selection_dictionary):
        '''Produce centroid arrays broken down by individual proteins for a given protein type and return in a dictionary categorized by protein type.'''
        dictionary_centroid_arrays = {}
        number_protein_atoms_dictionary = {'HA':3621,'NA':4094,'M2':372}
        for selection_descriptor, MDA_selection_object in current_MDA_selection_dictionary.iteritems():
            for protein_name in number_protein_atoms_dictionary.keys():
                if protein_name in selection_descriptor:
                    num_atoms_individual_protein_this_type = number_protein_atoms_dictionary[protein_name]
                    current_protein_name = protein_name
                    try:
                        current_index = 0
                        list_per_protein_centroids = []
                        while 1:
                            list_per_protein_centroids.append(MDA_selection_object[current_index: current_index + num_atoms_individual_protein_this_type].centroid())
                            current_index += num_atoms_individual_protein_this_type
                    except IndexError: #trying to lazily catch control flow when I've dealt with all protein atoms of a given type
                        pass
                    dictionary_centroid_arrays[current_protein_name] = numpy.array(list_per_protein_centroids)
        #print 'checkpoint 3: function dictionary keys:',dictionary_centroid_arrays.keys()
        return dictionary_centroid_arrays

    dict_MSD_values_this_replicate = {'MSD_value_dict':{},'MSD_std_dict':{},'frame_skip_value_list':[]} #for overall storage of MSD average / standard deviation values for this replicate

    for MSD_window_size_frames in MSD_window_size_list:
        list_per_window_average_displacements = []
        counter = 0
        trajectory_striding_dictionary = {} #store values in a dict as you stride through trajectory
        for ts in universe_object.trajectory[::MSD_window_size_frames]:
            #print 'checkpoint 4: MDA_full_protein_group_selection_dictionary.keys():',MDA_full_protein_group_selection_dictionary.keys()
            if counter == 0: #first parsed frame
                previous_frame_protein_centroid_array_dictionary = centroid_array_production(MDA_full_protein_group_selection_dictionary)
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            else: #all subsequent frames
                current_frame_protein_centroid_array_dictionary = centroid_array_production(MDA_full_protein_group_selection_dictionary)

                
                for protein_name in current_frame_protein_centroid_array_dictionary.keys():
                    if not protein_name in trajectory_striding_dictionary.keys(): 
                        trajectory_striding_dictionary[protein_name] = {'MSD_value_list_centroids':[]}
                        current_delta_array_centroids = previous_frame_protein_centroid_array_dictionary[protein_name] - current_frame_protein_centroid_array_dictionary[protein_name]
                        square_delta_array_centroids = numpy.square(current_delta_array_centroids)
                        sum_squares_delta_array_centroids = numpy.sum(square_delta_array_centroids,axis=1)
                        trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'].append(numpy.average(sum_squares_delta_array_centroids))
                    else:
                        current_delta_array_centroids = previous_frame_protein_centroid_array_dictionary[protein_name] - current_frame_protein_centroid_array_dictionary[protein_name]
                        square_delta_array_centroids = numpy.square(current_delta_array_centroids)
                        sum_squares_delta_array_centroids = numpy.sum(square_delta_array_centroids,axis=1)
                        trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'].append(numpy.average(sum_squares_delta_array_centroids))
                    #print 'checkpoint 1, trajectory_striding_dictionary.keys():',trajectory_striding_dictionary.keys()

                #reset the value of the 'previous' array as you go along:
                previous_frame_protein_centroid_array_dictionary = current_frame_protein_centroid_array_dictionary
                print multiprocessing.current_process().name, 'frame:', ts.frame 
            counter += 1
        #print 'checkpoint 2, trajectory_striding_dictionary.keys():',trajectory_striding_dictionary.keys()
        for protein_name,MSD_data_subdictionary in trajectory_striding_dictionary.iteritems():
            if not protein_name in dict_MSD_values_this_replicate['MSD_value_dict'].keys(): #initialize protein species subdictionaries as needed
                dict_MSD_values_this_replicate['MSD_value_dict'][protein_name] = []
                dict_MSD_values_this_replicate['MSD_std_dict'][protein_name] = []
                dict_MSD_values_this_replicate['MSD_value_dict'][protein_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_std_dict'][protein_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'])))
            else:
                dict_MSD_values_this_replicate['MSD_value_dict'][protein_name].append(numpy.average(numpy.array(trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'])))
                dict_MSD_values_this_replicate['MSD_std_dict'][protein_name].append(numpy.std(numpy.array(trajectory_striding_dictionary[protein_name]['MSD_value_list_centroids'])))
        dict_MSD_values_this_replicate['frame_skip_value_list'].append(MSD_window_size_frames)
    queue_IPC.put(dict_MSD_values_this_replicate)

def capture_plot_MSD_proteins(Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no'):
    '''Capture (and eventually plot) the per-protein MSD data.'''
    if from_pickle == 'no': #capture the data directly from the queue objects:
        sim35_MSD_tracking_dictionary_protein = Q2.get()
        sim36_MSD_tracking_dictionary_protein = Q3.get()
        sim37_MSD_tracking_dictionary_protein = Q4.get()
        sim38_MSD_tracking_dictionary_protein = Q5.get()
        sim39_MSD_tracking_dictionary_protein = Q6.get()
        sim_Danny_MSD_tracking_dictionary_protein = Q7.get()
        #pickle for safety / future plotting
        pickle.dump(sim35_MSD_tracking_dictionary_protein,open('sim35_MSD_protein.p','wb'))
        pickle.dump(sim36_MSD_tracking_dictionary_protein,open('sim36_MSD_protein.p','wb'))
        pickle.dump(sim37_MSD_tracking_dictionary_protein,open('sim37_MSD_protein.p','wb'))
        pickle.dump(sim38_MSD_tracking_dictionary_protein,open('sim38_MSD_protein.p','wb'))
        pickle.dump(sim39_MSD_tracking_dictionary_protein,open('sim39_MSD_protein.p','wb'))
        pickle.dump(sim_Danny_MSD_tracking_dictionary_protein,open('sim_Danny_MSD_protein.p','wb'))
    else: #load the data from pickled storage
        sim35_MSD_tracking_dictionary_protein,sim36_MSD_tracking_dictionary_protein,sim37_MSD_tracking_dictionary_protein,sim38_MSD_tracking_dictionary_protein,sim39_MSD_tracking_dictionary_protein,sim_Danny_MSD_tracking_dictionary_protein = [pickle.load(open(x,'rb')) for x in ['sim35_MSD_protein.p','sim36_MSD_protein.p','sim37_MSD_protein.p','sim38_MSD_protein.p','sim39_MSD_protein.p','sim_Danny_MSD_protein.p']]
        def register_diffusion_constant_plot_best_fit_line(output_dictionary, x_data_array, y_data_array, axis,index_window_filter):
            '''index_window_filter is for excluding values at short MSD time windows which may cause anomalous diffusion.'''
            x_data_array = x_data_array[index_window_filter:]
            y_data_array = y_data_array[index_window_filter:]
            slope, intercept = numpy.polyfit(x_data_array,y_data_array,1) #assuming that these are actually lines and not higher order polynomial trends
            #divide the slope (A ^ 2 / ns ) by 10 ^ 7 to convert to standard cm ^ 2 / s units:
            diffusion_constant = slope / (10 ** 7)
            #I also want to estimate the error in the diffusion constant (slope) values, and I'm going to mimic the approach of g_msd, which uses the difference between the slopes obtained for the two halves for the data:
            first_half_x_data, second_half_x_data = numpy.array_split(x_data_array,2)
            first_half_y_data, second_half_y_data = numpy.array_split(y_data_array,2)
            slope_first_half, intercept_first_half = numpy.polyfit(first_half_x_data,first_half_y_data,1) 
            slope_second_half, intercept_second_half = numpy.polyfit(second_half_x_data,second_half_y_data,1) 
            diffusion_constant_error_estimate = abs(slope_first_half - slope_second_half) / (10 ** 7)
            #put in the appropriate dictionary:
            output_dictionary[protein_label_string] = {'diffusion_constant':diffusion_constant,'error_estimate':diffusion_constant_error_estimate}
            #plot the best linear fit to the data for this particular lipid species:
            coefficients = slope, intercept
            polynomial = numpy.poly1d(coefficients)
            x = numpy.arange(0,600,100)
            ys = polynomial(x)
            axis.plot(x,ys,color='black') #best fit line

        def calculate_scaling_exponent_plot_best_fit_line(log_x_data_array,log_y_data_array,axis,protein_label):
            '''For calculation of scaling exponent (alpha) and plotting of the best fit line on the log plots.'''
            slope, intercept = numpy.polyfit(log_x_data_array,log_y_data_array,1) #assuming that these are actually lines and not higher order polynomial trends
            scaling_exponent = slope
            coefficients = slope, intercept
            polynomial = numpy.poly1d(coefficients)
            x = numpy.arange(0,3,0.1)
            ys = polynomial(x)
            axis.plot(x,ys,color='black') #best fit line
            #axis.text(0.5,0.5,'$\\alpha = {scaling_exponent}$ {protein_label}'.format(scaling_exponent = scaling_exponent, protein_label = protein_label))
            scaling_exponent_string = '%.2f' % scaling_exponent
            alpha_string = '$\\alpha = {scaling_exponent}$ {protein_label}'.format(scaling_exponent = scaling_exponent_string, protein_label = protein_label)
            return alpha_string





        def diffusion_constant_bar_chart(diffusion_constant_dictionary,title,axis):
            x_positions = numpy.arange(0,len(diffusion_constant_dictionary))
            color_list = ['blue','red','purple','black','green','cyan']
            bars_plotted = 0
            for lipid_label,diffusion_dictionary in diffusion_constant_dictionary.iteritems():
                #diffusion_constant = '%.3e' % diffusion_constant #convert to scientific notation
                axis.bar(x_positions[bars_plotted],diffusion_dictionary['diffusion_constant'],label=lipid_label,color=color_list[bars_plotted],width=0.7,yerr=diffusion_dictionary['error_estimate'],ecolor='black')
                print lipid_label, diffusion_dictionary['diffusion_constant']
                axis.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axis.legend(bbox_to_anchor=(0.5,-0.05),loc='upper center',numpoints=1,prop={'size':6},ncol=6,fancybox = True)
                bars_plotted += 1
            axis.set_title(title)
            axis.set_ylabel('D (cm$^2$ / s)')
            axis.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
            axis.set_xlim(0,6)
        import matplotlib.pyplot
        colordict = {'NA':'yellow', 'HA':'orange','M2':'pink'}
        fig = matplotlib.pyplot.figure()
        fig_log = matplotlib.pyplot.figure()

        subplot_number = 1
        for MSD_tracking_dictionary,plot_name in zip([sim35_MSD_tracking_dictionary_protein,sim36_MSD_tracking_dictionary_protein,sim37_MSD_tracking_dictionary_protein,sim38_MSD_tracking_dictionary_protein,sim39_MSD_tracking_dictionary_protein,sim_Danny_MSD_tracking_dictionary_protein],['virion 295 K', 'protein-restrained virion 295 K','virion 323 K','FORS virion 295 K','FORS virion 323 K','Danny viron (40% CHOL) 310 K']):
            centroid_diffusion_constant_dictionary = {}
            #ax = fig.add_subplot(6,2,subplot_number)
            ax_log = fig_log.add_subplot(6,1,subplot_number)
            print 'plot name:', plot_name
            #print 'MSD_tracking_dictionary.keys():',MSD_tracking_dictionary.keys()
            #print 'MSD_tracking_dictionary[MSD_value_dict].keys():', MSD_tracking_dictionary['MSD_value_dict'].keys()
            #print 'MSD_tracking_dictionary[MSD_std_dict].keys():', MSD_tracking_dictionary['MSD_std_dict'].keys()
            frame_window_array = numpy.array(MSD_tracking_dictionary['frame_skip_value_list'])
            if frame_window_array.max() == 500: #Danny's sim already skips 10 frames so convert directly to ns
                time_window_array_nanoseconds = frame_window_array
            else: #need to divide by 10 for other replicates
                time_window_array_nanoseconds = frame_window_array / 10.0
            for protein_label_string, MSD_centroid_value_array in MSD_tracking_dictionary['MSD_value_dict'].iteritems():
                #print 'protein_label_string:',protein_label_string
                #print 'MSD_centroid_value_array:', MSD_centroid_value_array
                #ax.errorbar(time_window_array_nanoseconds,MSD_centroid_value_array,yerr=MSD_tracking_dictionary['MSD_std_dict'][protein_label_string],label='{protein_label} centroids'.format(protein_label=protein_label_string),ls='none',marker='o')
                alpha_string = calculate_scaling_exponent_plot_best_fit_line(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_centroid_value_array),ax_log,protein_label_string)
                ax_log.scatter(numpy.log10(time_window_array_nanoseconds),numpy.log10(MSD_centroid_value_array),label='{alpha_string}'.format(alpha_string=alpha_string),marker='o',color=colordict[protein_label_string])
                #register_diffusion_constant_plot_best_fit_line(output_dictionary=centroid_diffusion_constant_dictionary, x_data_array=time_window_array_nanoseconds, y_data_array=MSD_centroid_value_array, axis=ax,index_window_filter=0)
                #ax.set_title(plot_name)
                ax_log.set_title(plot_name)
                #ax.set_ylim(0,3000)
                ax_log.set_ylim(0,4.0)
                #ax.set_xlim(0,550)
                ax_log.set_xlim(-0.1,3.0)
                #ax.legend(loc=2,numpoints=1,prop={'size':10})
                #ax.set_xlabel('Time (ns)')
                ax_log.set_xlabel('log(t)')
                ax_log.set_ylabel('log(MSD)')
            ax_log.legend(loc=2,numpoints=1,prop={'size':10})
            subplot_number += 1
            #ax2 = fig.add_subplot(6,2,subplot_number + 1)
            #diffusion_constant_bar_chart(diffusion_constant_dictionary=centroid_diffusion_constant_dictionary,title='D comparison',axis=ax2)
            #subplot_number += 2
        #fig.set_size_inches(14,29)
        #fig.savefig('MSD_tracking_flu_proteins.png',dpi=300)
        fig_log.set_size_inches(6,29)
        fig_log.savefig('log_scaling_exponent_MSD_tracking_flu_proteins.png',dpi=300)

def radial_distribution_lipid_tails_to_water(queue_IPC,universe_object):
    '''Influenza g(r) code for lipid tails to water.'''
    #slightly tricky implementation attempt: give 4 cores to this function so that the 7X dispatches on different cores will expand to use 28/32 cores
    #so this will actually become a 'sub-parent' of the module parent
    import collections
    def log_result_pool(result):
        parent_deque_rdfs.append(result)
        
    def per_core_work(lipid_tail_coordinates,solvent_coordinate_sub_array,bins,dmin,dmax):
        distance_matrix = scipy.spatial.distance.cdist(lipid_tail_coordinates,solvent_coordinate_sub_array)
        rdf_this_core, edges_this_core = numpy.histogram(distance_matrix,bins=nbins,range=(dmin,dmax))
        return rdf_this_core

    if universe_object.selectAtoms('resname DUPC').numberOfAtoms() > 0: #Danny's sim
        lipid_tail_selection_string = '(resname DUPC and name C4B) or (resname DPPC and name C4B)'
        frame_skip_value = 10 #this replicate already skips 10
    else:
        lipid_tail_selection_string = '(resname POPS and name C5B) or (resname DOPX and name C5B) or (resname DOPE and name C5B) or (resname FORS and name C4B) or (resname PPCH and name C4B)' #attempt to select the 'innermost' tail particles of each lipid species (exclude CHOL for fairer comparison with dengue system)
        frame_skip_value = 100

    solvent_selection_string = '(resname W or resname WF or resname ION) and around 100.0 ({lipid_tail_string})'.format(lipid_tail_string = lipid_tail_selection_string) #because flu is so massive, limit selected waters to those within a 100.0 A shell of the lipid tails

    lipid_tail_selection = universe_object.selectAtoms(lipid_tail_selection_string)
    solvent_selection = universe_object.selectAtoms(solvent_selection_string)
    n = lipid_tail_selection.numberOfAtoms() * solvent_selection.numberOfAtoms()

    dmin, dmax = 0.0, 110.0
    nbins = 200
    rdf, edges = numpy.histogram([0], bins=nbins, range=(dmin, dmax))
    rdf *= 0
    rdf = rdf.astype(numpy.float64)
    cumulative_boxvolume = 0
    frame_counter = 0
    for ts in universe_object.trajectory[::frame_skip_value]: #use a frame_skip_value to match the dengue skip 100 filtering
        pool = multiprocessing.Pool(4) #trying to work with a 'sub-Pool'; we'll see if this actually works
        parent_deque_rdfs = collections.deque() #place rdfs from individual cores here
        box_dimensions = universe_object.dimensions
        cumulative_boxvolume += MDAnalysis.coordinates.core.box_volume(box_dimensions)
        lipid_tail_coordinates = lipid_tail_selection.coordinates()
        solvent_coordinates = solvent_selection.coordinates()
        #try to split the array of solvent coordinates into 3000 sub-arrays that can be passed to different cores (trying to use cdist on smaller matrices over many cores vs. a single huge cdist job)
        list_solvent_coordinate_subarrays = numpy.array_split(solvent_coordinates,3000)
        for solvent_coordinate_sub_array in list_solvent_coordinate_subarrays:
            pool.apply_async(per_core_work,args=(lipid_tail_coordinates,solvent_coordinate_sub_array,nbins ,dmin,dmax),callback = log_result_pool)
        pool.close()
        pool.join()
        for sub_rdf in parent_deque_rdfs: #add in all the rdfs from individual cores in this frame
            rdf += sub_rdf
        frame_counter += 1
        print simulation_name, 'frame:', frame_counter

    numframes = frame_counter
    average_volume = cumulative_boxvolume / float(numframes)
    #normalize RDF
    radii = 0.5*(edges[1:] + edges[:-1])
    number_density = float(n) / average_volume
    vol = (4./3.)*numpy.pi*number_density*(numpy.power(edges[1:],3)-numpy.power(edges[:-1], 3))
    rdf = rdf / (vol * numframes)
    queue_IPC.put([radii,rdf]) #should be easy to determine which dataset belongs with which replicate based on the Q number


def capture_plot_radial_distribution(Q1,Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no'):
    '''For now, just pickle the captured data for safety / future plotting use.'''
    sim33_radial_dist_data = Q1.get()
    sim35_radial_dist_data = Q2.get()
    sim36_radial_dist_data = Q3.get()
    sim37_radial_dist_data = Q4.get()
    sim38_radial_dist_data = Q5.get()
    sim39_radial_dist_data = Q6.get()
    sim_Danny_radial_dist_data = Q7.get()
    #pickle for safety / future plotting
    pickle.dump(sim33_radial_dist_data,open('sim33_radial_dist_data.p','wb'))
    pickle.dump(sim35_radial_dist_data,open('sim35_radial_dist_data.p','wb'))
    pickle.dump(sim36_radial_dist_data,open('sim36_radial_dist_data.p','wb'))
    pickle.dump(sim37_radial_dist_data,open('sim37_radial_dist_data.p','wb'))
    pickle.dump(sim38_radial_dist_data,open('sim38_radial_dist_data.p','wb'))
    pickle.dump(sim39_radial_dist_data,open('sim39_radial_dist_data.p','wb'))
    pickle.dump(sim_Danny_radial_dist_data,open('sim_Danny_radial_dist_data.p','wb'))

def radial_distribution_protein_types_to_lipids(queue_IPC,universe_object,testing=None):
    '''Calculate the RDFs between protein TMD centroids and surrounding lipids.'''
    import scipy.spatial.distance
    dictionary_data_current_replicate = {}
    #use lipid centroids
    lipid_selection_dict = {'POPS':{} , 'DOPX':{}, 'DOPE':{}, 'FORS':{}, 'PPCH': {} ,'DUPC': {},'DPPC': {},'CHOL':{} } #should cover my sims + Danny's
    for lipid_name, sub_dictionary in lipid_selection_dict.iteritems():
        lipid_selection_string = 'resname {res_string}'.format(res_string = lipid_name)
        selection = universe_object.selectAtoms(lipid_selection_string)
        lipid_selection_dict[lipid_name]['selection_string'] = lipid_selection_string
        lipid_selection_dict[lipid_name]['MDA_selection_object'] = selection
    lipid_dict_copy = lipid_selection_dict.copy()
    for lipid_name, sub_dictionary in lipid_selection_dict.iteritems():
        selection = lipid_selection_dict[lipid_name]['MDA_selection_object'] 
        if selection.numberOfResidues() == 0:
            del lipid_dict_copy[lipid_name] #to avoid exceptions for empty groups in distance matrices (i.e., Danny's vs. my sims have different lipids from the set above)
    lipid_selection_dict = lipid_dict_copy #reassign to the filtered dict
    lipid_selection_dict_final_100ns = lipid_selection_dict.copy() #I think trying to recycle the lipid_selection_dict for the final 100 ns analysis may be causing issues in the final 100 ns g(r) calculations--it's at least worth trying the use of a fresh copy for the second iteration of the g(r) analysis
    debug_dict_3 = lipid_selection_dict.copy() #debug
    #I'll have to pull out the residue objects and their centroids within the trajectory iteration for loop
    #now, for the proteins the selection is a bit tricky because I may want to select only the TMD regions as the ectodomains extend quite far from the bilayer proper; also, how am I going to treat the palmitoyl chains on HA? There is some relevant code in protein_proximity_analysis_v2() above, which can probably be recycled here.

    HA_selection = universe_object.selectAtoms('bynum 1:289680') #80 HA molecules * 3621 particles / molecule
    NA_selection = universe_object.selectAtoms('bynum 289681:338808') #12 NA molecules * 4094 particles / molecule
    M2_selection = universe_object.selectAtoms('bynum 338809:344388') #15 M2 molecules * 372 particles / molecule

    def TMD_particle_selector(input_array,molecule_type):
        '''Selects the TMD coordinate elements from the input_array and combines to a simplified new numpy array with TMD particle coordinates only. The input_array should be a single numpy array of coordinates for a single protein of molecule_type.'''
        if molecule_type == 'HA': #the index numbers are based on study of topology combined with Danny's DPhil thesis
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

    def initialize_histogram_data_lipid_subdictionaries(lipid_subdictionary,dmin,dmax,nbins):
        for protein_string in ['HA','NA','M2']:
            rdf, edges = numpy.histogram([0], bins=nbins, range=(dmin, dmax))
            rdf *= 0
            rdf = rdf.astype(numpy.float64)
            lipid_subdictionary['{protein_string}_rdf'.format(protein_string=protein_string)], lipid_subdictionary['edges'] = rdf, edges

    def traj_parser_rdf_protein_lipid(start_frame,end_frame,iteration_label,input_dict,lipid_dict):
        dmin, dmax = 0.0, 25.0 #25.0 A shell limit
        nbins = 100
        cumulative_boxvolume = 0
        frame_counter = 0
        record_testing_data = 'no'
        for ts in universe_object.trajectory[start_frame:end_frame:10]: #every 10th frame = 1 ns intervals
            print multiprocessing.current_process().name, 'frame:', ts.frame, 'label:', iteration_label
            box_dimensions = universe_object.dimensions
            cumulative_boxvolume += MDAnalysis.coordinates.core.box_volume(box_dimensions)
            HA_coordinates,NA_coordinates,M2_coordinates = [selection.coordinates() for selection in [HA_selection,NA_selection,M2_selection]]
            list_individual_HA_protein_coordinate_arrays = numpy.split(HA_coordinates,80) #split to list of coord arrays for each of the 80 HA molecules
            list_individual_NA_protein_coordinate_arrays = numpy.split(NA_coordinates,12) #similarly for NA
            list_individual_M2_protein_coordinate_arrays = numpy.split(M2_coordinates,15) #similarly for M2
            assert list_individual_HA_protein_coordinate_arrays[0].shape == (3621,3), "Incorrect shape for HA coordinates in list_individual_HA_protein_coordinate_arrays."
            assert list_individual_NA_protein_coordinate_arrays[0].shape == (4094,3), "Incorrect shape for NA coordinates in list_individual_NA_protein_coordinate_arrays."
            assert list_individual_M2_protein_coordinate_arrays[0].shape == (372,3), "Incorrect shape for M2 coordinates in list_individual_M2_protein_coordinate_arrays."
            list_HA_TMD_coordinate_arrays = [TMD_particle_selector(HA_coord_array,'HA') for HA_coord_array in list_individual_HA_protein_coordinate_arrays]
            list_NA_TMD_coordinate_arrays = [TMD_particle_selector(NA_coord_array,'NA') for NA_coord_array in list_individual_NA_protein_coordinate_arrays]
            list_M2_TMD_coordinate_arrays = [TMD_particle_selector(M2_coord_array,'M2') for M2_coord_array in list_individual_M2_protein_coordinate_arrays]
            array_HA_TMD_centroids = numpy.array([numpy.average(HA_TMD_array,axis=0) for HA_TMD_array in list_HA_TMD_coordinate_arrays])
            array_NA_TMD_centroids = numpy.array([numpy.average(NA_TMD_array,axis=0) for NA_TMD_array in list_NA_TMD_coordinate_arrays])
            array_M2_TMD_centroids = numpy.array([numpy.average(M2_TMD_array,axis=0) for M2_TMD_array in list_M2_TMD_coordinate_arrays])
            assert array_NA_TMD_centroids.shape == (12,3), "Incorrect shape for array_NA_TMD_centroids: %s" % array_NA_TMD_centroids.shape
            if record_testing_data == 'no': #store the coordinates of all protein TMD centroids and all phosphate particles for later inspection (to confirm that I'm getting reasonable protein TMD centroids relative to the lipid 'halo' later on when I'm doing validation / analysis of the pickled results in IPython notebook)
                input_dict['all_protein_TMD_centroids_array_first_frame_check_{label}'.format(label=iteration_label)] = numpy.concatenate((array_HA_TMD_centroids,array_NA_TMD_centroids,array_M2_TMD_centroids))
                input_dict['phosphate_halo_array_first_frame_check_{label}'.format(label=iteration_label)] = universe_object.selectAtoms('name PO4').coordinates()
                record_testing_data = 'yes' #only record from the first frame of current iteration
            #now for the lipids:
            for lipid_name, subdictionary in lipid_dict.iteritems():
                #print subdictionary.keys()
                if not subdictionary.has_key('HA_rdf'): #initialize the empty rdf histogram data if it doesn't exist yet
                    #print 'initializing new RDF histograms' #**debug problem--this is only getting called 3 times on the first frame of first traj_parser function call 
                    initialize_histogram_data_lipid_subdictionaries(lipid_subdictionary = subdictionary,dmin=dmin,dmax=dmax,nbins=nbins)
                current_lipid_residue_object_list = lipid_dict[lipid_name]['MDA_selection_object'].residues
                current_lipid_array_residue_centroids = numpy.array([residue.centroid() for residue in current_lipid_residue_object_list])
                for subdict_rdf_key, protein_TMD_centroid_array in zip(['HA_rdf','NA_rdf','M2_rdf'],[array_HA_TMD_centroids,array_NA_TMD_centroids,array_M2_TMD_centroids]):
                    #print 'current_lipid_array_residue_centroids.shape',current_lipid_array_residue_centroids.shape #(11487,) -- so that's a problem in current debug --it's the correct number of DPPC centroids (Danny's sim) but it's a flat array
                    #print 'protein_TMD_centroid_array.shape', protein_TMD_centroid_array.shape
                    distance_matrix = scipy.spatial.distance.cdist(current_lipid_array_residue_centroids,protein_TMD_centroid_array)
                    current_rdf, current_edges = numpy.histogram(distance_matrix,bins=nbins,range=(dmin,dmax))
                    #if lipid_name == 'CHOL': #selective debugging
                        #print iteration_label, lipid_name, subdict_rdf_key, 'histogram max current frame:', current_rdf.max() #debugging
                    #subdictionary[subdict_rdf_key] += current_rdf
                    subdictionary[subdict_rdf_key] = subdictionary[subdict_rdf_key] + current_rdf
                    #if lipid_name == 'CHOL': #selective debugging
                        #print iteration_label, lipid_name, subdict_rdf_key, 'current summed RDF max:', subdictionary[subdict_rdf_key].max() #debugging
                    subdictionary['n_value_' + subdict_rdf_key] = current_lipid_array_residue_centroids.shape[0] * protein_TMD_centroid_array.shape[0]   #the 'n' parameter used in rdf analysis / normalization is done separately for each lipid / protein pair
            frame_counter += 1
        #after trajectory iteration, deal with the normalization / final calculations for each lipid-protein rdf data set in the dictionary
        numframes = frame_counter
        average_volume = cumulative_boxvolume / float(numframes)
        edges = numpy.histogram([0], bins=nbins, range=(dmin, dmax))[1] #can probably avoid this dummy histrogram just to pull out the edges, but ok for now
        #now make the iterative changes in the dictionary for normalization of rdf:
        for lipid_name, subdictionary in lipid_dict.iteritems():
            #print lipid_name
            subdictionary['numframes'] = numframes #debug
            subdictionary['cumulative_boxvolume'] = cumulative_boxvolume #debug
            subdictionary['average_volume'] = average_volume #debug
            for subdict_rdf_key in ['HA_rdf','NA_rdf','M2_rdf']:
                number_density = float(subdictionary['n_value_' + subdict_rdf_key]) / average_volume
                vol = (4./3.)*numpy.pi*number_density*(numpy.power(edges[1:],3)-numpy.power(edges[:-1], 3))
                subdictionary[subdict_rdf_key] = subdictionary[subdict_rdf_key] / (vol * numframes) #should be the appropriately normalized RDF data for this particular lipid:protein pair
                #if lipid_name == 'CHOL': #selective debugging for output sanity
                    #print iteration_label, lipid_name, subdict_rdf_key, 'number density:', number_density, 'vol.shape:', vol.shape, 'vol.max():', vol.max(),'max rdf:', subdictionary[subdict_rdf_key].max(),'non-zero min rdf:', subdictionary[subdict_rdf_key][subdictionary[subdict_rdf_key] > 0].min()
        for lipid_name, subdictionary in lipid_dict.iteritems():
            input_dict[iteration_label + lipid_name] = {}
            for key, value in subdictionary.iteritems():
                if not key == 'MDA_selection_object': #these selection objects can't be pickled so I don't want them moving between processes in queues (it will fail)
                    input_dict[iteration_label + lipid_name][key] = value #should eventually get dumped back to queue
        for lipid_name, subdictionary in lipid_dict.iteritems():
            for key in subdictionary.keys():
                if key in ['HA_rdf','NA_rdf','M2_rdf','n_value_M2_rdf','edges','n_value_HA_rdf','n_value_NA_rdf']:
                    del subdictionary[key] #I seem to be having problems with carry-over of these RDF dictionary entries between calls to traj_parser_rdf_protein_lipid--trying to deal with that

    traj_parser_rdf_protein_lipid(start_frame=1,end_frame=1000,iteration_label='first_100_ns',input_dict=dictionary_data_current_replicate,lipid_dict=lipid_selection_dict)
    traj_parser_rdf_protein_lipid(start_frame=-1000,end_frame=-1,iteration_label='final_100_ns',input_dict=dictionary_data_current_replicate,lipid_dict=lipid_selection_dict_final_100ns)
    

    if not testing:
        queue_IPC.put(dictionary_data_current_replicate) 
    else:
        return dictionary_data_current_replicate


def capture_RDF_data_protein_types_to_lipids(Q2,Q3,Q4,Q5,Q6,Q7,from_pickle = 'no'):
    '''Just capture / pickle the RDF data--I plan to dissect the pickled data in IPython.'''
    sim35_protein_lipid_RDF_data = Q2.get()
    sim36_protein_lipid_RDF_data = Q3.get()
    sim37_protein_lipid_RDF_data = Q4.get()
    sim38_protein_lipid_RDF_data = Q5.get()
    sim39_protein_lipid_RDF_data = Q6.get()
    sim_Danny_protein_lipid_RDF_data = Q7.get()
    for sim_name,RDF_data in zip(['sim35','sim36','sim37','sim38','sim39','sim_Danny'],[sim35_protein_lipid_RDF_data,sim36_protein_lipid_RDF_data,sim37_protein_lipid_RDF_data,sim38_protein_lipid_RDF_data,sim39_protein_lipid_RDF_data,sim_Danny_protein_lipid_RDF_data]):
        pickle.dump(RDF_data,open('{sim_name}_protein_lipid_RDF_data.p'.format(sim_name = sim_name),'wb'))


def capture_RDF_data_protein_types_to_lipids_testing(Q7,from_pickle = 'no'): #temporarily only accept Q7 for small test purposes
    sim_Danny_protein_lipid_RDF_data = Q7.get()
    pickle.dump(sim_Danny_protein_lipid_RDF_data,open('{sim_name}_protein_lipid_RDF_data.p'.format(sim_name = 'sim_Danny'),'wb'))


def code_to_dispatch(queue_object,coordinate_file_path,ordered_list_of_trajectory_paths,analysis_function):
    '''Code to execute on individual cores (child processes). Will always create a Universe object first and then perform the analysis on this Universe object.Will optionally use the queue_object to communicate data back to the parent process if there is something that can be pickled that should be communicated (Universe object and file descriptors can't easily be communicated this way). Simply change the analysis function to perform different analyses on the trajectories in parallel.'''
    print multiprocessing.current_process().name, 'Starting' #to monitor start of this process on terminal
    #always generate the MDA Universe object first:
    universe_object = generate_merged_MDA_universe(coordinate_file_path,ordered_list_of_trajectory_paths)
    #now execute the analysis function for this universe object on this core:
#    analysis_function(universe_object) #the form of this f'n may change depending on number of arguments accepted, etc.
    analysis_function(queue_object, universe_object) #the form of this f'n may change depending on number of arguments accepted, etc.
    #analysis_function(queue_object, universe_object,'PPCH') #the form of this f'n may change depending on number of arguments accepted, etc.
    print multiprocessing.current_process().name, 'Finishing' #to monitor end of this process on terminal





def main():
    '''Main control function of code.'''
    start_time = time.time() #the start time in seconds so I can benchmark

    #sims 35/36 (see Table in my lab notes) share the same coordinate file for MDA topology purposes:
    #sim 37 should also share this coordinate file as it is simply the equivalent of sim35 (295K) at 323K
    sim_35_36_coordinate_file_compact_no_solvent = '/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/July2_2012_equil_vesicle_NO_Forss/protein_incorporation/g_membed/completed_embedding_Aug14_2012/compact_sys_no_solvent.gro'

    #create properly ordered lists (part 1, part 2, part 3, part ...) of trajectory files for each of the sims I want to parse:
    #I'm now starting to make adjustments such that more recent data stored on /sansom/n04/bioc1009/ can also be included
    list_sim33_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/July2_2012_equil_vesicle_NO_Forss/','*no_solvent*xtc'),key= lambda file_string: int(file_string[99:101].replace('_',''))) #sort by file name part number
    list_sim33_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories('/sansom/n04/bioc1009/July2_2012_equil_vesicle_NO_Forss/','*no_solvent*xtc'),key= lambda file_string: int(file_string[-30:-28])))
    list_sim35_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/Aug_20_2012_non_Forss_virion_prod_5us_rad_080_unrestrained/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_',''))) #idea with the latter slicing is to leave room for part numbers with two digits to be properly sorted
    list_sim36_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/Aug_21_2012_non_Forss_virion_prod_5us_rad_080_PULL_CODE/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_','')))
    list_sim36_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories('/sansom/n04/bioc1009/Aug_21_2012_non_Forss_virion_prod_5us_rad_080_PULL_CODE/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_',''))))
    list_sim37_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/Nov_26_2012_non_Forss_virion_prod_5us_rad_080_unrestrained_323K/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_','')))
    list_sim38_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/Jan_11_2013_sim38_FORS_inclusive_295K_unrestrained/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_','')))
    list_sim38_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories('/sansom/n04/bioc1009/Jan_11_2013_sim38_FORS_inclusive_295K_unrestrained/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_',''))))
    list_sim39_trajectories_compact_no_solvent = sorted(produce_list_trajectories('/sbcb/hfs0/bioc1009/postdoc_data/production_simulations/March_16_2013_sim39_FORS_inclusive_323K_unrestrained/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_','')))
    list_sim39_trajectories_compact_no_solvent.extend(sorted(produce_list_trajectories('/sansom/n04/bioc1009/March_16_2013_sim39_FORS_inclusive_323K_unrestrained/','*no_solvent*xtc'),key = lambda file_string: int(file_string[-6:-4].replace('_',''))))
    list_sim_Danny_trajectories_compact_no_solvent = ['/sbcb/hfs0/parton/flu/virion_build/a-5000ns/parts0-39-noW-skip10.xtc'] #just using a single trajectory file produced by Danny, but note the SKIP 10 processing!

    #test proper ordering of trajectory lists when needed:
    #trajectory_list_order_checking(list_sim_Danny_trajectories_compact_no_solvent) #for whichever list you need to check
    #import sys; sys.exit('debug exit for trajectory order checking')

    #reduce trajectory sizes when doing debug/testing work:
    testing = 'no'
    #testing = 'yes'
    if testing == 'yes':
        list_sim33_trajectories_compact_no_solvent = list_sim33_trajectories_compact_no_solvent[0:1]
        list_sim35_trajectories_compact_no_solvent = list_sim35_trajectories_compact_no_solvent[0:1]
        list_sim36_trajectories_compact_no_solvent = list_sim36_trajectories_compact_no_solvent[0:1]
        list_sim37_trajectories_compact_no_solvent = list_sim37_trajectories_compact_no_solvent[0:1]
        list_sim38_trajectories_compact_no_solvent = list_sim38_trajectories_compact_no_solvent[0:1]
        list_sim39_trajectories_compact_no_solvent = list_sim39_trajectories_compact_no_solvent[0:1]

    #sim35 has 52 .xtc files for loading and is hanging into a zombie process when I try to load it, so I am going to do some debugging with systematic variation of the number of .xtc files allowed for sim35:
    #list_sim35_trajectories_compact_no_solvent = list_sim35_trajectories_compact_no_solvent[0:40]


    #now I'm basically going to open MDA Universe objects and perform analyses on these objects on separate cores and I will avoid interprocess communication (i.e., from children back to this parent) until I have something useful to communicate that can be pickled (this is a requirement--so can't simply communicate the MDA Universe objects back to the parents)

    jobs = [] #for tracking parallel processes

    #generate a set of queue objects as these will be the first arguments to the code_to_dispatch() function:
    queue_1, queue_2, queue_3, queue_4, queue_5,queue_6,queue_7 = (multiprocessing.Queue(), multiprocessing.Queue(),multiprocessing.Queue(),multiprocessing.Queue(), multiprocessing.Queue(),multiprocessing.Queue(),multiprocessing.Queue())
    
    #----change this to perform a different parallel analysis---:
    #current_parallel_analysis = frame_iterator
    #current_parallel_analysis = track_PPCH_PO4_OD
    #current_parallel_analysis = track_multiple_protein_COGs_in_virion
    #current_parallel_analysis = leaflet_asymmetry_analysis
    #current_parallel_analysis = leaflet_asymmetry_analysis_2
    #current_parallel_analysis = stratification_tracking
    #current_parallel_analysis = positional_probability_virus_species
    #current_parallel_analysis = protein_proximity_analysis
    #current_parallel_analysis = protein_proximity_analysis_v2
    #current_parallel_analysis = lipid_protein_closest_approach_analysis
    #current_parallel_analysis = euclidean_displacement_assessment
    #current_parallel_analysis = sphericity_tracking
    #current_parallel_analysis = mean_square_displacement
    #current_parallel_analysis = mean_square_displacement_by_species
    #current_parallel_analysis = mean_square_displacement_proteins
    #current_parallel_analysis = radial_distribution_lipid_tails_to_water 
    current_parallel_analysis =  radial_distribution_protein_types_to_lipids

    #now, each of the simulations will be opened and parsed on a separate core via a multiprocessing.Process object with some arguments passed to the code_to_dispatch(queue_object,coordinate_file_path,ordered_list_of_trajectory_paths,analysis_function) function:
    sim_33_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_1,'/sansom/sc2/bioc1009/postdoc_work/improved_virion/vesicle_build_adjusted_June19_2012/compact_neutral_antifreeze_nosolvent.gro',list_sim33_trajectories_compact_no_solvent,current_parallel_analysis))
    sim_35_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_2,sim_35_36_coordinate_file_compact_no_solvent,list_sim35_trajectories_compact_no_solvent,current_parallel_analysis))
    sim_36_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_3,sim_35_36_coordinate_file_compact_no_solvent,list_sim36_trajectories_compact_no_solvent,current_parallel_analysis))
    sim_37_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_4,sim_35_36_coordinate_file_compact_no_solvent,list_sim37_trajectories_compact_no_solvent,current_parallel_analysis))
    sim_38_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_5,'/sansom/sc2/bioc1009/postdoc_work/improved_virion/Forss_virion_by_conversion_Sept_2012/md_Lizzie_PR3_compact.gro',list_sim38_trajectories_compact_no_solvent,current_parallel_analysis))
    sim_39_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_6,'/sansom/sc2/bioc1009/postdoc_work/improved_virion/Forss_virion_by_conversion_Sept_2012/md_Lizzie_PR3_compact.gro',list_sim39_trajectories_compact_no_solvent,current_parallel_analysis)) #can use same coord file as sim38 b/c sim39 is the same coords at an elevated temperature
    sim_Danny_process = multiprocessing.Process(target = code_to_dispatch, args = (queue_7,'/sbcb/hfs0/parton/flu/virion_build/a-5000ns/prod-noW.tmp.gro',list_sim_Danny_trajectories_compact_no_solvent,current_parallel_analysis)) #for Danny's 1.8 microsecond virion simulation using 40% CHOL along with DPPC & DUPC; note the SKIP 10 filtering of this trajectory relative to my own trajectories

    #comment these blocks to turn off the parallel MDA parsing and just go straight to plotting--------------------
    #now start the multiprocessing.Process objects on separate cores:
    #for job in [sim_33_process,sim_35_process,sim_36_process,sim_37_process,sim_38_process,sim_39_process,sim_Danny_process][1:]: 
        #jobs.append(job)
        #job.start() #start loading each of the MDA Universes into memory from different cores and then perform the specified analysis code

    #now use a data capture and plotting function to deal with the parallel output of the child processes (usually matches the analysis function):
    #capture_plot_PPCH_PO4_OD(queue_1, queue_2, queue_3,queue_4,queue_5,queue_6,from_pickle = 'yes')
    #capture_plot_multiple_protein_COGs_in_virion(queue_2,queue_3,from_pickle = 'yes')
    #capture_plot_leaflet_asymmetry_analysis(queue_1,queue_2,queue_3,from_pickle = 'yes')
    #capture_plot_leaflet_asymmetry_analysis_2(queue_1,queue_2,queue_3,from_pickle = 'yes')
    #capture_plot_stratification_tracking(queue_1,queue_2,queue_3,queue_4,queue_5,queue_6,from_pickle = 'no')
    #capture_plot_positional_probability_virus_species(queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'yes')
    #capture_plot_protein_proximity_analysis(queue_2,queue_3,queue_4,queue_5,from_pickle = 'yes')
    #capture_plot_protein_proximity_analysis_v2(queue_2,queue_3,queue_4,queue_5,from_pickle = 'yes')
    #capture_plot_lipid_protein_closest_approach_analysis(queue_2,queue_3,queue_4,queue_5,from_pickle = 'yes')
    #capture_plot_euclidean_displacement_assessment(queue_1,queue_2,queue_3,queue_4,queue_5,from_pickle = 'yes',lipid_name = 'PPCH')
    capture_plot_sphericity_tracking(queue_1,queue_2,queue_3,queue_4,queue_5,queue_6,from_pickle = 'yes')
    #capture_plot_MSD(queue_1,queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'yes',sim_danny_only='no')
    #capture_plot_MSD_by_lipid(queue_1,queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'yes')
    #capture_plot_MSD_proteins(queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'yes')
    #capture_plot_radial_distribution(queue_1,queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'no')
    #plot_log_MSD_by_lipid_scaling_exponent()
    #capture_RDF_data_protein_types_to_lipids(queue_2,queue_3,queue_4,queue_5,queue_6,queue_7,from_pickle = 'no')
    #capture_RDF_data_protein_types_to_lipids_testing(queue_7,from_pickle = 'no')

    #now block control flow in this (parent) script until the children processes have completed:
    #so, this code will probably work best (relative to serial analysis) for cases where each of the trajectories is parsed for a similar time
    #note that the join() methods will cause a deadlock if placed ahead of the above capture functions that use the .get() method of Queue
    #for job in jobs:
        #job.join() #join method will block indefinitely until the child process is done
    #----------------------------------------------------------------------------------------------------------------------
    print 'All children have exited' #debugging

    end_time = time.time() #the end time in seconds
    execution_time_seconds = end_time - start_time
    execution_time_minutes = execution_time_seconds / 60.0
    execution_time_hours = execution_time_minutes / 60.0
    print 'execution time (s): ', execution_time_seconds, 'execution time (minutes):', execution_time_minutes, 'execution time (hours):', execution_time_hours

if __name__ == '__main__': #execute the main control function only if this file is called as a top-level script
    main()
