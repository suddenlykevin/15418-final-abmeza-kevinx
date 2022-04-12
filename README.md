
# 15-418 Final Project Proposal: Parallel Pixelated Image Abstraction
Anthony Meza (abmeza), Kevin Xie (kevinx), 15-418 (S22)

# Summary
We are going to parallelize [Gerstner et al.'s approach](https://gfx.cs.princeton.edu/pubs/Gerstner_2012_PIA/Gerstner_2012_PIA_full.pdf) to the pixelization/downsampling and color quantization of an image while retaining its salient details. The main focus will be recreating and parallelizing each iteration of the image processing algorithm and yielding an image which resembles pixel art.

# Background

Pixelated image abstraction is a term coined in Gerstner et al.'s research paper of the same name, used to describe the process of abstracting high resolution images into low resolution outputs that resemble pixel art. In essence, the term refers to a method of creating an output image that is comparable to that of the image manually recreated by a pixel artist.

The two primary constraints for this process are the limited size of the color palette and the lower resolution of the output, both of which are user defined. In order to address these two constraints, Gerstner’s algorithm uses Simple Linear Interactive Clustering (SLIC) and Mass Constrained Deterministic Annealing (MCDA). Both are clustering algorithms which together are able to segment the “super pixels” in the image and decide on the colors within the pallet, respectively. 

The general algorithm used will be the following, as proposed by Gerstner et al.:

 - initialize superpixels, palette and temperature T (Section 4.1) .
 - while (T > Tf ) .
   - refine superpixels with 1 step of modified SLIC
   - associate superpixels to colors in the palette
   - refine colors in the palette
   - if (palette converged)
     - reduce temperature T = α T
     - expand palette
 - post-process

Below is a visualization of the process, which involves setting up the image, and then afterwards a repetition of (c) assigning pixels to “super pixels”, which are associated with the pixels in our new image, and (d) assigning superpixels to the colors in our palette.


# The Challenge

Graphics and image filtering is a classically parallelizable task. Reading and rendering the I/O files should be easily parallelizable. However, compared to some more naive approaches, this filtering algorithm iterates over multiple distinct phases that could be a challenge to parallelize effectively. 

Performing SLIC, for example, requires associating a cluster of neighboring pixels to a specific superpixel, which changes upon every iteration. This will require figuring out how to group out pixels for each super pixel, and managing groups of various different sizes for consecutive iterations. Following each segmentation, Gerstner et al. also smooths the segment boundaries using Laplacian smoothing, which will require substantial communication between superpixels.

It will also likely be a challenge to try and find ways to optimize efficient memory accessing algorithms given the variability in the shape and size of the clusters for each iteration and how they could be generalized over variable threads. 

There might also be potential for communication between threads that might need to be considered. For example, in the expanding portion of the algorithm, clusters are split into different color blocks, which could be more efficiently implemented with communication between different threads to decide which pixels are associated with which group.

# Resources
We intend to write the software in C++ according to the specifications set forth in Gerstner et al.’s original paper. As the algorithm looks to be graphics-intensive and quite parallelizable, we plan on using relying mainly on CUDA to make the program parallel.

In addition to the original paper, Gerstner et al.’s paper heavily cites [Achanta et al.’s SLIC method for superpixel segmentation](https://infoscience.epfl.ch/record/149300/files/SLIC_Superpixels_TR_2.pdf) and [Rose et al.’s Mass Constrained Deterministic Annealing](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3047&rep=rep1&type=pdf) (which we may need to simplify) for clustering and palette refinement.

# Goals and Deliverables

## Plan to Achieve
 - A sequential implementation which matches or closely approximates Gerstner et al.’s original filtering behavior
   - This should be completed early on to leave enough time for experimenting with parallelization strategies.
 - Parallel superpixel segmentation
   - We believe this could be a good dimension for parallelism because it 
 - Refine palette in parallel over superpixels
   - We believe this could be a good dimension for parallelism because the operation is data parallel between superpixels. Each palette refinement operation takes place over a single superpixel segment.
 - Render final image in parallel
   - We believe this would be the easiest parallelism to implement using a similar approach to renderCircles in Assignment 2.
 - Brief analysis and comparison

## Hope to Achieve
 - User-defined palette color selections
   - In addition to limiting the number of colors that the palette can contain, it may be interesting to find a way to constrain exactly which colors the filter can use. This would forgo the MCDA palette refinement in favor of dynamic color warping/quantization [as described by Kim et al. here](https://dl.acm.org/doi/10.1145/3450626.3459776).
 - User defined importance maps (as described in Gerstner et al.)
   - Gerstner et al. briefly describe a strategy for users to define the relative importance of regions of an image using an “importance map”, then using the map to dedicate the majority of the color palette’s dynamic range to the more “important” regions of the image.

## Stretch
 - Extend to video processing
   - It could be cool to extend past static images to processing the frames of a video in order to convert an entire video to resemble pixel-art animation.
 - Adapt as photoshop plug-in
   - Since Photoshop plug-ins are written against a C++ API, this could be a feasible way to allow others to more easily access and experiment with our solution

# Platform Choice
We have chosen to use C++, leveraging CUDA on the RTX 2080 GPUs of the **GHC machines**. Because the GHC machine’s GPU has the most CUDA cores of Carnegie Mellon’s local resources, we think it would be the best way to develop and test our implementation. We would also like to eventually test our implementation on the **Bridges-2 supercomputer** to test how performance scales with CUDA cores available.

# Planned Schedule

**Week 1** - Read papers, scout resources online, and begin writing sequential implementation.<br>
- [x] Read papers
- [x] Scout resources

**Week 2** - Complete sequential implementation and focus on parallelizing superpixel segmentation.<br>
- [x] Week 1 overflow: Begin writing sequential implementation
- [ ] Complete sequential implementation
- [ ] Focus on parallelizing superpixel implementation

**Week 3** - Complete superpixel segmentation, begin to parallelize palette refinement and image render.<br>
- [ ] Complete superpixel parallelization
- [ ] Begin to parallelize palette refinement and image render

**Week 4** - Complete parallelization, test on Bridges-2 and begin writing the final report.<br>
- [ ] Complete parallelization
- [ ] Test on bridges-2
- [ ] Begin final report

**Week 5** - Complete and submit the final report, prepare for the presentation.<br>
- [ ] Complete and submit final report
- [ ] Prepare presentation

# Milestone Report
In our planned schedule, we stated that by the second week, we should have completed our sequential implementation and begun some preliminary experimentation with parallelization. Unfortunately, it seems we were a little bit too optimistic about our comprehension of the underlying algorithm of pixelation image abstraction. As a result, we decided to spend the first week focused on better understanding the algorithms cited in the Gerstner paper. From this, we gained a better understanding of SPIC superpixel segmentation and MCDA. It is clear that the main source of complexity for our project, as well as the dimensions of parallelization, will come from the multiple stages that are involved with the algorithm we are trying to build. 

With respect to our goals and deliverables, we are currently a little behind schedule, since we are still trying to finish the sequential implementation of our project. It has proven to be a little more challenging since there are a lot of new concepts and formulas that we have been forced to analyze in order to figure out exactly how to proceed. We do not think we will be able to accomplish our stretch goals, but we remain confident that we should still be able to achieve our minimum deliverables at the end of our project. Since we have been able to develop a solid understanding of our algorithm, progress from here should be manageable.

Our focus at the moment is managing the rounding of values which has been an issue that we have been dealing with when converting between color spaces. This is especially important since color data in our algorithm requires a decent level of precision. We have also been looking into taking advantage of image libraries to speed up the prototyping process for the sequential implementation. This is, of course, a minor speed bump as we dive into implementing the full algorithm.

## Updated Schedule
**Week 3**<br>
- [ ] Complete sequential implementation
- [ ] Parallelize simple tasks (colorspace conversion, image read/render)
- [ ] Think about superpixel parallelization

**Week 4**<br>
- [ ] Benchmark sequential implementation (Bridges-2)
- [ ] Parallelize superpixel segmentation
- [ ] Begin parallelizing palette refinement and laplacian smoothing
- [ ] Begin final report

**Week 5**<br>
- [ ] Complete parallelization
- [ ] Benchmark parallel implementation (Bridges-2)
- [ ] Complete and submit final report
- [ ] Prepare presentation

## Poster session
At our poster session we plan on being able to show a demo of our working project running both sequentially and in parallel to demonstrate the speedup we were able to achieve, and in the process demonstrate what our code is able to accomplish. In addition to this, we will hopefully be able to show some diagrams with statistics of the speed up across a different number of processors. 
