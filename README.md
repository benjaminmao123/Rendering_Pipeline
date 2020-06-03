# Rendering Pipeline

A basic rendering pipeline. 

The pipeline implements 5 steps of a typical rendering pipeline.

1. Implements a vertex shader which processes each vertex into an output vertex.

2. Assembles triangle primitives.

3. Performs clipping, perspective divide, and viewport transformation to window space.

4. Interpolation of each vertex which are broken down into fragments to be processed by a fragment shader.

5. Depth test using z-buffer for hidden surfaces.

6. Rasterization 
