# Hyperbolic Orbifold Tutte Embeddings

Matlab code implementing the [Siggraph Asia 2016 paper, "Hyperbolic Orbifold Tutte Embeddings"](https://noamaig.github.io/html/projects/hyperbolic/hyperbolic_low.pdf).

An extension of Tutte's embedding to hyperbolic domains, namely hyperbolic orbifolds which enable seamless parameterization and mapping.


### Example scripts:
- `example_embed_mesh_sphere` - map a spherical mesh to a hyperbolic orbifold.
- `example_embed_mesh_disk` - map a disk mesh to a hyperbolic orbifold.
- `example_map` - map meshes to one another by embedding both into the same hyperbolic orbifold.

### Main classes:
- `flattener` - the object used to compute the embedding. 
- `mapper` - object that takes a bunch of flatteners and computes the mapping induced by them.

Please note the optimization scheme used in this implementation is more advanced than in the original paper, as it uses L-BFGS optimization with Laplacian preconditioning in the spirit of ["Accelerated quadratic proxy for geometric optimization"](https://dl.acm.org/citation.cfm?id=2925920). For additional details, see appendix B in ["Spherical orbifold Tutte embeddings"](https://dl.acm.org/citation.cfm?id=3073615) wherein the same implementation was used.

The code is provided as-is for academic use only and without any guarantees. Please contact the author to report any bugs.
Written by [Noam Aigerman](https://noamaig.github.io/) and [Shahar Kovalsky](https://github.com/shaharkov).

 
