import matplotlib
matplotlib.use('svg')
from pylab import *
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import math
import matplotlib.pyplot as plt
from ReadDepth import ReadDepth
from mRNAsObject import mRNAsObject


def plot_density_single(read_depth_object, sample_label,
                        mRNAs, strand,
                        graphcoords, graphToGene, axvar,
                        min_circRNA_depth,min_exon_junction,
                        sum_circRNA_expression,
                        paired_end=False,
                        intron_scale=30,
                        exon_scale=4,
                        color='r',
                        ymax=None,
                        number_junctions=True,
                        # resolution=.5,
                        showXaxis=True,
                        showYaxis=True,
                        nyticks=3,
                        nxticks=4,
                        show_ylabel=True,
                        show_xlabel=True,
                        font_size=6,
                        numbering_font_size=6,
                        junction_log_base=10,
                        plot_title=None,
                        plot_label=None):

    def junc_comp_function(a,b):

        '''
            junc_comp_function is a __cmp__ function which allows junctions to be
                sorted based on the length of the intron that they span

                junctions with shorter intronic regions come first
        '''
        
        a_coordinates = map(int, a.split(':')[1].split('-'))
        b_coordinates = map(int, b.split(':')[1].split('-'))

        a_distance = a_coordinates[1] - a_coordinates[0]
        b_distance = b_coordinates[1] - b_coordinates[0]

        return a_distance - b_distance
    
    # extract data from read_depth_object
    tx_start = read_depth_object.low
    tx_end = read_depth_object.high
    chrom = read_depth_object.chrm
    wiggle = read_depth_object.wiggle
    jxns = read_depth_object.junctions_dict
    circRNA_ = read_depth_object.circRNA
    # print circRNA_
    # print sum_circRNA_expression
    if sum_circRNA_expression == False:
        circRNA = sorted(circRNA_.items(), key=lambda item:item[1])
    else:
        if circRNA_ == {}:
            circRNA = []
        else:
            circRNA_sorted = sorted(circRNA_.items(), key=lambda item:item[1])
            # print circRNA_sorted
            name = circRNA_sorted[0][0].split(':')[0]
            circRNA_all = []
            for idx in range(len(circRNA_sorted)):
                temp_circRNA = [0 for i in range(4)]
                temp_circRNA[0] = int(circRNA_sorted[idx][0].split(':')[1].split('-')[0])
                temp_circRNA[1] = int(circRNA_sorted[idx][0].split(':')[1].split('-')[1])
                temp_circRNA[2] = circRNA_sorted[idx][1]
                temp_circRNA[3] = -1
                circRNA_all.append(temp_circRNA)
            for i in range(len(circRNA_all)):
                start_idx = circRNA_all[i][0]
                end_idx = circRNA_all[i][1]
                if circRNA_all[i][3] != -1:
                    continue
                circRNA_all[i][3] = i
                for j in range(i+1, len(circRNA_all)):
                    if circRNA_all[j][0] == start_idx and circRNA_all[j][3] == -1:
                        circRNA_all[j][3] = i
                for k in range(i+1, len(circRNA_all)):
                    if circRNA_all[k][1] == end_idx and circRNA_all[k][3] == -1:
                        circRNA_all[k][3] = i
            # print circRNA_all
            for flag in range(len(circRNA_all)):
                value_sum = 0.0
                for idx in range(len(circRNA_all)):
                    if circRNA_all[idx][3] == flag:
                        value_sum += circRNA_all[idx][2]
                for idxx in range(len(circRNA_all)):
                    if  circRNA_all[idxx][3] == flag:
                        circRNA_all[idxx][2] = value_sum
            # print circRNA_all
            circRNA = []
            index = 0
            for idx_ in range(len(circRNA_all)):
                temp = []
                tempName = name + ':' + str(circRNA_all[idx_][0]) +'-'+ str(circRNA_all[idx_][1])
                temp.append(tempName)
                temp.append(circRNA_all[idx_][2])
                circRNA.append(temp)
                index += 1

    print circRNA

    maxheight = max(wiggle)
    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.5 * ymax   

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    prevx = graphcoords[0]
    tmpval = []
    for i in range(len(graphcoords)):
        tmpval.append(wiggle[i])
        # if abs(graphcoords[i] - prevx) > resolution:
        compressed_wiggle.append(mean(tmpval))
        compressed_x.append(prevx)
        prevx = graphcoords[i]
        tmpval = []

    fill_between(compressed_x, compressed_wiggle,\
        y2=0, color=color, lw=0)
   
    sslists = []
    for mRNA in mRNAs:
        tmp = []
        for s, e in mRNA:
            tmp.extend([s, e])
        sslists.append(tmp)

    # sort the junctions by intron length for better plotting look
    jxns_sorted_list = sorted(jxns.keys(),cmp=junc_comp_function)
    current_height = -3 * ymin / 4
    # plot junction on the top part
    for plotted_count, jxn in enumerate(jxns_sorted_list):
        # print plotted_count, jxn, jxns[jxn]
        leftss, rightss = map(int, jxn.split(":")[1].split("-"))

        ss1, ss2 = [graphcoords[leftss - tx_start - 1],\
            graphcoords[rightss - tx_start]]

        mid = (ss1 + ss2) / 2

        label_x = (ss1+ss2)/2

        leftdens = wiggle[leftss - tx_start - 1]
        rightdens = wiggle[rightss - tx_start]

        pts = [(ss1, leftdens),
                ((ss1+ss2)/2, (leftdens +rightdens + 2*current_height)/2),
                (ss2, rightdens)]
        midpt = cubic_bezier(pts, .5)
        label_y = (leftdens +rightdens + 2*current_height)/2

        # print min_exon_junction
        if number_junctions and jxns[jxn] >= min_exon_junction:
            text(label_x, label_y, '{0}'.format(round(jxns[jxn],2)),
                 fontsize=numbering_font_size, ha='center', va='center', backgroundcolor='w')
        
        # a = Path(pts, [Path.MOVETO, Path.LINETO, Path.LINETO])
        # set minimum junction_log_base
        if min_exon_junction >= 1.0 and jxns[jxn] >= min_exon_junction:
            a = Path(pts, [Path.MOVETO, Path.LINETO, Path.LINETO])
            p = PathPatch(a, ec=color, lw=log(jxns[jxn] + 1) / log(junction_log_base), fc='none')
            axvar.add_patch(p)
        elif min_exon_junction < 1.0:
            if jxns[jxn] > min_exon_junction and jxns[jxn] <= 1.0:
                temp_junction_log_base = 1.3
                a = Path(pts, [Path.MOVETO, Path.LINETO, Path.LINETO])
                p = PathPatch(a, ec=color, lw=log(jxns[jxn] + 1) / log(temp_junction_log_base), fc='none')
                axvar.add_patch(p)
            elif jxns[jxn] > 1.0:
                a = Path(pts, [Path.MOVETO, Path.LINETO, Path.LINETO])
                p = PathPatch(a, ec=color, lw=log(jxns[jxn] + 1) / log(junction_log_base), fc='none')
                axvar.add_patch(p)

    # plot circRNA on the bottom part
    for plotted_count, rna in enumerate(circRNA):
        # print plotted_count, rna
        leftss, rightss = map(int, rna[0].split(":")[1].split("-"))
        ss1, ss2 = [graphcoords[leftss - tx_start - 1], graphcoords[rightss - tx_start]]

        mid = (ss1 + ss2) / 2
        label_x = (ss1+ss2)/2

        pts = [(ss1, 0), (ss1, -current_height*1.2), (ss2, -current_height*1.2), (ss2, 0)]
        midpt = cubic_bezier4(pts, .5)
        label_y = -current_height

        if number_junctions and rna[1] >= min_circRNA_depth:
            text(label_x, label_y, '{0}'.format(round(rna[1],2)),
                 fontsize=numbering_font_size, ha='center', va='center', backgroundcolor='w')
        # print min_circRNA_depth
        if min_circRNA_depth >= 1.0 and rna[1] >= min_circRNA_depth:
            a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            p = PathPatch(a, ec=color, lw=log(rna[1] + 1) / log(junction_log_base), fc='none')
            axvar.add_patch(p)
        elif min_circRNA_depth < 1.0:
            if rna[1] > min_circRNA_depth and rna[1] <= 1.0:
                temp_junction_log_base = 1.3
                a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                p = PathPatch(a, ec=color, lw=log(rna[1] + 1) / log(temp_junction_log_base), fc='none')
                axvar.add_patch(p)
            elif rna[1] > 1.0:
                a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                p = PathPatch(a, ec=color, lw=log(rna[1] + 1) / log(junction_log_base), fc='none')
                axvar.add_patch(p)

    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')

    if showXaxis:
        axvar.xaxis.set_ticks_position('bottom')
        xlabel('Genomic coordinate (%s), "%s" strand'%(chrom, strand), fontsize=font_size)

        max_graphcoords = max(graphcoords) - 1
        xticks(linspace(0, max_graphcoords, nxticks),
               [graphToGene[int(x)] for x in \
                linspace(0, max_graphcoords, nxticks)],
               fontsize=font_size)
    else:
        axvar.spines['bottom'].set_color('none')
        xticks([])

    xlim(0, max(graphcoords))
    return axvar


# Plot density for a series of bam files.
def plot_density(settings,event,read_depths_dict,mRNA_object,ordered_genotypes_list,genotype_dict_counts):

    intron_scale = settings["intron_scale"]
    exon_scale = settings["exon_scale"]
    colors = settings["colors"]
    ymax = settings["ymax"]
    number_junctions = settings["number_junctions"]
    # resolution = settings["resolution"]
    junction_log_base = settings["junction_log_base"]
    reverse_minus = settings["reverse_minus"]
    font_size = settings["font_size"]
    nyticks = settings["nyticks"]
    nxticks = settings["nxticks"]
    show_ylabel = settings["show_ylabel"]
    show_xlabel = settings["show_xlabel"]
    plot_title = settings["plot_title"]
    numbering_font_size = settings["numbering_font_size"]
    min_circRNA_depth = settings["min_circRNA_depth"]
    min_exon_junction = settings["min_exon_junction"]
    sum_circRNA_expression = settings["sum_circRNA_expression"]
    # print min_circRNA_depth,min_exon_junction
    # Always show y-axis for read densities for now
    showYaxis = True
    
    # parse mRNA_object to get strand, exon_starts, exon_ends, tx_start, tx_end, chrom
    strand = mRNA_object.strand
    chom = mRNA_object.chrm
    exon_starts = mRNA_object.exon_starts
    exon_ends = mRNA_object.exon_ends
    tx_start = mRNA_object.low
    tx_end = mRNA_object.high
    mRNAs = mRNA_object.mRNAs

    # Get the right scalings
    graphcoords, graphToGene = getScaling(tx_start, tx_end, strand,
                                          exon_starts, exon_ends, intron_scale,
                                          exon_scale, reverse_minus)
    # print graphcoords

    nfiles = len(read_depths_dict.keys())

    if plot_title is not None and plot_title != '':
        # Use custom title if given
        suptitle(plot_title, fontsize=10)
    elif plot_title == '':
        suptitle(event, fontsize=10)
        
    plotted_axes = []

    labels_list = []

    for i, group_genotype in enumerate(ordered_genotypes_list):
        average_read_depth = read_depths_dict[group_genotype]

        if colors is not None:
            color = colors[i]
        else:
            color = None
        if i < nfiles - 1:
            showXaxis = False 
        else:
            showXaxis = True 

        ax1 = subplot2grid((nfiles + 2, 1), (i,0), colspan=1)
        
        # Read sample label
        sample_label = group_genotype
        labels_list.append(group_genotype+' (n = '+str(genotype_dict_counts[group_genotype])+')')

        plotted_ax = plot_density_single(read_depth_object=average_read_depth,
                        sample_label=sample_label,
                        mRNAs=mRNAs, strand=strand,
                        graphcoords=graphcoords,graphToGene=graphToGene,axvar=ax1,
                        paired_end=False,
                        intron_scale=intron_scale,
                        exon_scale=exon_scale,
                        color=color,
                        ymax=ymax,
                        number_junctions=number_junctions,
                        # resolution=resolution,
                        showXaxis=showXaxis,
                        showYaxis=showYaxis,
                        nyticks=nyticks,
                        nxticks=nxticks,
                        show_ylabel=show_ylabel,
                        show_xlabel=show_xlabel,
                        font_size=font_size,
                        numbering_font_size=numbering_font_size,
                        junction_log_base=junction_log_base,
                        min_circRNA_depth=min_circRNA_depth,
                        min_exon_junction=min_exon_junction,
                        sum_circRNA_expression=sum_circRNA_expression)

        plotted_axes.append(plotted_ax)

    ##
    ## Figure out correct y-axis values
    ##
    ymax_vals = []
    if ymax != None:
        # Use user-given ymax values if provided
        max_used_yval = ymax
    else:
        # Compute best ymax value for all samples: take
        # maximum y across all.
        used_yvals = [curr_ax.get_ylim()[1] for curr_ax in plotted_axes]
        # Round up
        max_used_yval = 10*math.ceil(max(used_yvals)*0.1)

    # Reset axes based on this.
    # Set fake ymin bound to allow lower junctions to be visible
    fake_ymin = -0.5 * max_used_yval
    universal_yticks = linspace(0, max_used_yval, nyticks + 1)

    # Round up yticks
    universal_ticks = map(math.ceil, universal_yticks)
    for sample_num, curr_ax in enumerate(plotted_axes):
        if showYaxis:
            curr_ax.set_ybound(lower=fake_ymin, upper=1.2*max_used_yval)
            curr_yticklabels = []
            for label in universal_yticks:
                if label <= 0:
                    # Exclude label for 0
                    curr_yticklabels.append("")
                else:
                    if label % 1 != 0:
                        curr_yticklabels.append("%.1f" %(label))
                    else:
                        curr_yticklabels.append("%d" %(label))
            curr_ax.set_yticklabels(curr_yticklabels, fontsize=font_size)
            curr_ax.spines["left"].set_bounds(0, max_used_yval)
            curr_ax.set_yticks(universal_yticks)
            curr_ax.yaxis.set_ticks_position('left')
            curr_ax.spines["right"].set_color('none')
            if show_ylabel:
                y_horz_alignment = 'left'
                curr_ax.set_ylabel('Mean depth', fontsize=font_size, va="center", ha=y_horz_alignment,labelpad=10)

        else:
            curr_ax.spines["left"].set_color('none')
            curr_ax.spines["right"].set_color('none')
            curr.ax.set_yticks([])
        ##
        ## Plot sample labels
        ##
        sample_color = colors[sample_num]
        # Make sample label y position be halfway between highest
        # and next to highest ytick
        if len(universal_yticks) >= 2:
            halfway_ypos = (universal_yticks[-1] - universal_yticks[-2]) / 2.
            label_ypos = universal_yticks[-2] + halfway_ypos
        else:
            label_ypos = universal_yticks[-1]
        curr_label = labels_list[sample_num]
        label_ypos = label_ypos * 1.2
        curr_ax.text(max(graphcoords), label_ypos, curr_label, 
            fontsize=font_size, va='bottom', ha='right', color=sample_color)
                

    # Draw gene structure
    ax1 = subplot2grid((nfiles + 2, 1), (nfiles,0), colspan=1,rowspan=2)
    plot_mRNAs(tx_start, mRNAs, strand, graphcoords, reverse_minus)
    subplots_adjust(hspace=.1, wspace=.7)


def getScaling(tx_start, tx_end, strand, exon_starts, exon_ends,
               intron_scale, exon_scale, reverse_minus):
    """
    Compute the scaling factor across various genic regions.
    """
    # print tx_start, tx_end, strand
    # print exon_starts, exon_ends
    # print intron_scale, exon_scale, reverse_minus

    exoncoords = zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start : exon_ends[i] - tx_start] = 1

    graphToGene = {}
    graphcoords = zeros((tx_end - tx_start + 1), dtype='f')
    x = 0
    if strand == '+' or not reverse_minus:
        for i in range(tx_end - tx_start + 1):
            graphcoords[i] = x
            graphToGene[int(x)] = i + tx_start
            if exoncoords[i] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    else:
        for i in range(tx_end - tx_start + 1):
            graphcoords[-(i + 1)] = x
            graphToGene[int(x)] = tx_end - i + 1
            if exoncoords[-(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    return graphcoords, graphToGene


def plot_mRNAs(tx_start, mRNAs, strand, graphcoords, reverse_minus):
    """
    Draw the gene structure.
    """
    yloc = 0 
    exonwidth = .3
    narrows = 50

    for mRNA in mRNAs:
        # print mRNA
        for s, e in mRNA:
            s = s - tx_start
            e = e - tx_start
            # print graphcoords[s], graphcoords[e]
            x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]
            y = [yloc - exonwidth / 2, yloc - exonwidth / 2,\
                yloc + exonwidth / 2, yloc + exonwidth / 2]
            fill(x, y, 'k', lw=.5, zorder=20)

        # Draw intron.
        #axhline(yloc, color='k', lw=.5)
        plot([min(graphcoords),max(graphcoords)],[yloc,yloc], color='k',lw=0.5)

        # Draw intron arrows.
        spread = .2 * max(graphcoords) / narrows
        for i in range(narrows):
            loc = float(i) * max(graphcoords) / narrows
            if strand == '+' or reverse_minus:
                x = [loc - spread, loc, loc - spread]
            else:
                x = [loc + spread, loc, loc + spread]
            y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
            plot(x, y, lw=.5, color='k')

        yloc += 1 

    xlim(0, max(graphcoords)) 
    ylim(-.5, len(mRNAs) + .5)
    box(on=False)
    xticks([])
    yticks([]) 


def cubic_bezier(pts, t):

    p0, p1, p2 =pts
    p0 = array(p0)
    p1 = array(p1)
    p2 = array(p2)
    return p0 * (1 - t)**3 + 3 * t * p1 * (1 - t) ** 2 + t**3 * p2

def cubic_bezier4(pts, t):
    """
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = array(p0)
    p1 = array(p1)
    p2 = array(p2)
    p3 = array(p3)
    return p0 * (1 - t)**3 + 3 * t * p1 * (1 - t) ** 2 + 3 * t**2 * (1 - t) * p2 + t**3 * p3

def draw_sashimi_plot(output_file_path,settings,var_pos,average_depths_dict,mRNAs_object,ordered_genotypes_list,genotype_dict_counts):

    '''
        draw_sashimi_plot draws the complete sashimi plot

        output_file_path is the file path that the plot will be written to

        settings is a dict containing the settings for the sashimi plot

        var_pos is the location of the SNP, in the format chr1:12345

        average_depths_dict is a dict containing the average read depths by genotype. The keys are the genotypes,
            and the values are ReadDepth objects

        mRNAs_object is an mRNAsObject containing information about the transcript structure

        plot_title is the title of the plot

        return values:
            None. Draws sashimi plot

    '''

    plt.figure(figsize=[settings['width'],settings['height']])
    plot_density(settings,var_pos,average_depths_dict,mRNAs_object,ordered_genotypes_list,genotype_dict_counts)
    plt.savefig(output_file_path,transparent=True)
