"""
  Low level code generator
  Uses the MTK backend.
Generates
1. An unoptimized easy to read representation of a MTK system
Or
2. A low level DifferentialEquations.jl based system for further processing or C code generation.

Code generated by this code generator are expected to be written to file.
See top level API for how to call this.
It should also be noted that this code generator is a bit more experimental than the main MTK one.
Also since it does not subdivide the equations, for large 10,000+ systems the compilation time might be higher.

TODO:
Integrate the external code generator
"""
module SystemsExporter


end