# Background Information
### Important Terms
* *__Cycle:__* each “layer” of the DNA sequencing (one chemical cycle of adding dNTPs (addition of a nucleotide to each DNA read), washing away excess nucleotides, and recording fluorescence from the most-recently-added nucleotides)  
* *__Tile:__* one predetermined area of the flow cell (think of it separated into a grid of “tiles”)  
  * When multiple cycles are run, one “tile” extends down in a 3D column through all cycles of that area of the grid, so collectively a “tile” represents this 3D column  
  * One tile contains many DNA reads (distinct clusters that have been detected by the sequencer)  
  * Currently the code assigns one tile of information to each reader to process (the readers run in parallel)  
* *__Lanes:__* a flow cell is separated into lanes of DNA sequencing that sometimes has multiple rows of tiles  
* *__Parts:__* sequencing data from one lane is stored in multiple “parts” labeled with the name of that lane due to the sheer size of the data  
  * An abstract term used only to account for the way the data is stored, no actual separation in the “parts” on a flow cell  
* *__Quality score:__* the sequencer’s confidence that the nucleotide it recorded at a given cluster was actually the nucleotide that was added  
  * Diversity of reads generally increases quality score: when there are many different nucleotides being sequenced during one cycle, there is no overwhelming fluorescence of one nucleotide that might cause interference of nearby clusters (reason for adding PhiX to sequencing sample)  
    * *__Diversity of reads:__* having an equal representation of each nucleotide at each cycle
  * Overclustering decreases quality score (more difficult to distinguish between nearby clusters if there is high saturation of the flow cell)  
  * Quality score generally decreases as length of read increases since more and more reads within each cluster get “out of sync” in their nucleotide addition and the clusters begin to get harder to distinguish  
* *__Filter (aka Cluster Pass Filter):__* after 25 cycles, the sequencer makes a decision of whether a specific read is faulty or not; if it has too many nearby clusters to distinguish or is just perceived as an empty area of the flow cell, it will no longer be tracked by the sequencer  
  * A filter is a binary mask that is generated after these initial cycles and applies to all future cycles
  * One filter for each tile or can consider combining all of the individual tile filters to create one massive flow cell filter

![](https://sfvideo.blob.core.windows.net/sitefinity/images/default-source/product-page-images/next-generation-sequencing/ngs_adapter_designs.png?sfvrsn=8ce20807_8)  
source and additional info: https://www.idtdna.com/pages/products/next-generation-sequencing/adapters  
* *__Index (aka barcode):__* usually two 8bp indices per sample, one for each end of a paired-end read (blue and orange portions of figure above)  
  * Separating the reads by these indices is essential at the end of data processing to separate the DNA reads of distinct experiments, that have been loaded onto the same flow cell, from each other  
* *__Adapter sequence:__* the actual sequence of a specific index (the “index” is just the label to identify it by)  
* *__UMIs:__* short sequences used to tag distinct molecules before PCR amplification and then identify each molecule after PCR (green portion of the figure above)  
  * The bcl2fastq2 conversion software removes UMIs from reads
* *__BCL:__* base call file  
* *__CBCL:__* our input files; concatenated base call files which contain aggregated data from BCL files of a run  
  * One cbcl file for each group of tiles from the same lane and surface
* *__Multiplexing:__* adding indices to each distinct sample/experiment so that they can be separated in post-sequencing analysis  
* *__Demultiplexing:__* assigning a read of a cluster to a sample/experiment based on its indices
  
### More on Gene Sequencing
* Video from Illumina on flow cell setup and bridge amplification (shows how each read is done and where each index ends up in the process): https://www.youtube.com/watch?v=fCd6B5HRaZ8&t=204s  
* Diagrams of Illumina’s NovaSeq read process:  
  * A thorough step-by-step:  
  ![](https://supportassets.illumina.com/content/dam/illumina-support/images/bulletins/PEcell1.png)  
  * Skips the bridge amplification steps but includes full read1 --> index1 --> index2 --> read2 pipeline: 
  ![](https://supportassets.illumina.com/content/dam/illumina-support/images/bulletins/PEcell2.png)  
* NovaSeq:  
  * Uses red/green wavelengths combos to assign a color to a base rather than 4 distinct base colors
