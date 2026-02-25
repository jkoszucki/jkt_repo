# Chapter2

Table: 'ktypes.csv' (TBD)
Table: 'ktypes_modifications.csv' (TBD)
Table: 'ktypes_sim.csv' (TBD)


# Plots API

#### Plot1: cumlative plot of ktypes
Cumulative plot of ktypes structures collected over the years. 
X: 'Year', Y: cumulative number of k-type structures collected

Split dataset based on 'Source': 
'KPAM' -> color:blue, label: k-pam databse
'LIT' -> color:red, label: literature


#### Plot2: modifications and monosaccharide residues prevalence

There are 2 subplots here which share the Y axis (xlim 0.05-1.05).
To plot this table you will need to load both tables with ktypes ('sample/ktypes.csv') and the table with modifications ('sample/modifications.csv'). Read their headers to learn the data structure before plotting.

Subplot 1. Shows the prevalence of 'pyruvylation' and 'acetylation' across the ktypes in a single stacked bar. 
Subplot 1 - 1st stacked bar:
Shows ratio of ktypes carrying 'pyruvylation' or 'acetylation' colored by bond type. When both modifications are present there are counter as category (Ac + Pyr). A layer correspoding to Ktypes without pyruvylation is 'gray'. 

The stacks order is sorted accordinlgy to the type of the bond a 'pyruvylation' is carrying, where the darkest color is on the bottom, then ligher colors, gray color (ktypes without pyruvylation). The same patter applys to 'acetylation', but in yellow. Where both modifications are present then as red. Importantly, the bond type is not always known, to these assign separate category i.e. 'unknown'

The order of increasing darnkess of the yellow color for 'acetylation': 0-6, O-2, O-, unknown
The order of increasing darnkess of the blue color for 'pyruvylation': O-2,3 -> O-3,4 -> O-4,3 -> O-4,6 -> O-6,4, unknown

Subplot 2. Shows the prevalence of monosaccharide residues in ktypes along with distribution of 'pyruvylation' and 'acetylation' across the monosaccharide residues.
Each column corresponds to unique monosaccharides residues and its prevalence in ktypes. If the monosaccharide carries the 'acetylation' or 'pyruvylation' there is a stack corresponding to the ratio of monosaccharides with this type of modification, relative to total number of modifications. There are separate stacks for each type of modification depending on the bond it carries. 

An example calculation:
total # of ktypes: 100
ktypes carrying Man residues: 50
Man residues carrying O-6-acetylation: 10
Man residues carrying O-2,4-pyruvylation: 20

Prevalence of the Man with O-6-acetylation: 10/100 -> stack 10% of thickness colored bright yellow
Prevalence of the Man with O-2,3-pyruvylation: 20/100 (20%) -> stack 20% of thickness colored brigh blue
Prevalence of the Man without modification: (50-10-20)/100 -> stack 20% of thickness colored gray 
Remaining 50% is white. For visibility draw outline of of the bar.


### Plot 3

Plot 3 shows combinations of all monosaccharides found in the ktypes along with their frequency.  and the subplot 2 shows the frequency of the given combination.

Load ktypes.csv table and from be 'backbone_linkages_coarse' and 'branch' extract all unique monosaccharides. Importantly, remove edges and suffix '-OUT' in 'backbone_linkages_coarse' and anchor from 'branch'.

For instance, in the 'backbone_linkages_coarse' cell GalA(α1-3)Man(α1-2)Man(α1-3)Gal(β1,2)-OUT one need to strip '-OUT' suffix and remove all edges in the brackets eg, (β1,2). In the 'branch' column Man(α1-4)GalA(P1) one would need to remove the acnhor i.e. 'GalA(P1)' and the edges in the brackets.

Please, implement this as a using .split and .strip methods as regex are difficult. Specifically, strip '-OUT' in backbone columns, then split by '(' and take every second element except the last one (edge). Then take the set to take unique values. After doing this for all ktypes draw the count plot sorted by frequency in the descending order.

These subplots share Y axis, where the left plot shows what are the unique monosaccharides along with their combinations and the right plot shows the frequency of each unique combination.
Subplot 1. Shows the combinations of unique monosaccharides as a presence/asbence matrix compositon. 
Subplot 2. Count plot of unique combinations of K-type compositions sorted by frequency.




# Draw API
Script: ktypes_draw.py.
Task: implementation of the modifications of the monosaccharide residues in the drawings.
Instruction:
Read information about the modifications use the ‚modifications’ table. 

The modifications are supposed to be undirected edges with „the ball” at the end. Specifically it supposed to be circle 1/3 of size of the regular circle radius along with the outline which is  half of the size of the regular line width. The edge of these modifications is half of the thickness of the regular edge. Accordingly to the type of the modification the color will vary: O-acetylation (OAc) will have red ball, O-pyruvylation (OPy) will have blue ball, and any other type of modification will have a gray ball. The label for the node is supposed to be above the node. Do not put the label for the edge.  

A dictionary to use between modification type and the label above the edge.
pyruvylation - Pyr
acetylation - Ac
formylation - For
lactylation - Lac
glutemate - Glu


Modifications frequency can vary. Accordingly to the frequency of the modification use different transparency of the ball (circle). Ranges of frequency to use:
<30 - very opaque
30-60 - slioghtly opaque
60-90 -  very slightly opaque
100 - dark
unknown fequency - gray


This codes the location of the modifications in the structure of the polysaccharide (CPS). Let’s consider few examples.  type=K1, modification=pyruvylation, monosaccharide=GlcA@C1
In the k-type K1 the modification pyruvylation is located at GlcA monosaccharide at position 1 of the core.


type=K2, modification=pyruvylation, monosaccharide=Man@B
In the k-type K1 the modification pyruvylation is located at Man monosaccharide in the branch. The monosaccharides have always only single monosaccharide of a given type, so there is no ambiguity where to put the modification. That’s why there is no position where there is branch with one exception (below).

The exception:
type=K2, modification=formylation, monosaccharide=Glc@B1
type=K2, modification=formylation, monosaccharide=Glc@B2

There are 2 branches in this case, each with single Glc, both of them contain formulation.

In the visualization I would like also include the modifications that don’t have known location i.e. in the modification column the cell has a value, but the monosaccharide cell empty. Could you please, then add, with the same style and coloring, an edge and node representing unassigned modification on the X: max(xlim)-0.5, max(ylim) + 0.5; add the 0.5 as to separate params to control the location of the free modification.

