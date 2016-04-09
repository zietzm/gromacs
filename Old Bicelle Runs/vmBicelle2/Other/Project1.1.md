# Lipid Bicelle Using Jiang Protocol Exactly

We have received a number of erroneous results. This can be for a number of reasons,
but in order to maintain the hoped-for simplicity of generating a lipid bicelle
by exactly following a methodology already developed, we are going to start over.

## Setup System

We first want to create a mixed-lipid bilayer. Jiang used a DPPC bilayer developed
by the Martini group directly.

Once we have the file, we need to replicate it laterally. This means that we
basically double the size and number of lipids inside.
```
genconf -f dppc_bilayer.gro -o double_bilayer.gro -nbox 2
```

This is a "lateral" replication, as Jiang puts it, meaning that the box is replicated
in the z-direction.
