# Bifurcation diagram of two populations of heterogeneous Kuramoto-Sakaguchi oscillators

These are commands for the Python interface of `auto-07p`.

At the end of this file we provide instructions to install the software, althought this might not work depending on your system configuration. See the official documentation.

Notice that $\mu$ in the manuscript corresponds to `p` in these codes.

## Provided files

The software requires at least two files to run, which we provide and should not be modified:

- `oa2.f90` provides the system equations and the corresponding Jacobian. We also modified the options in the file to provide the maxima and minima of limit-cycles. 
- `c.oa2` provides the default constants for `auto-07`. These can be overwritten within the program, thus no need to modify this file.  

## Bifurcations along $K$ for fixed $p=0.9$

### Fixed points:

First, we use Euler method to find a fixed point for $K=7$:

```
init=run('oa2',ICP=[1,2,3],IPS=-2,NMX=100000,PAR={'K': 7.0, 'p' : 0.9, 'alpha' : 1.2})
```
The fixed point corresponds to an assymetric (chimera) state, as $R_a\neqR_b$. 
We can continue this solution increasing $K$:

```
ic=init(201)
asym=run(ic,IPS=1,NMX=10000,ICP=[1,2])
```
We see that the continuation stops at $K=200$, because we specified so in the `c.oa2` file.
Also, a Hopf bifurcation (`HB`) has been detected at $K\approx 7.37$. We will investigate this later.
For the moment, we save it on a new variable:

```
hopf = asym('HB1')
```

The Python interface of `auto-07p` is just Python. Thus we can plot the previous results using `matplotlib`:

```
import matplotlib.pyplot as plt
plt.plot(asym['K'],asym['Ra'],asym['K'],asym['Rb'])
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
``` 
*(Depending on your version of numpy, this might produce errors due to deprecation of `np.bool`.
If this happens, a quick workaround is to do `import numpy as np; np.bool=np.bool_;`)*
We can now run backwards in $K$ by specifying `DS="-"`:

```
new_branch = run(ic,IPS=1,NMX=10000,ICP=[1,2],DS="-")
```

A bifurcation at $K\approx 6.66$ is detected (and then Auto turns backwards again). 
From this bifurcation, a new solution branch is detected, and auto computes the solution directly.

The two branches can be accessed as:
```
asym = new_branch[0] # The chimera state for the full $K$-range.
sym  = new_branch[1] # This corresponds to the homogeneous solution of the system.
```

Let's visualize the results:

```
plt.plot(asym['K'],asym['Ra'],'black')
plt.plot(asym['K'],asym['Rb'],'black')
plt.plot(sym['K'],sym['Ra'],'black')
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
```

We see that the homogeneous branch has been continued to negative values of $R$!
In order to obtain a physically meaningfull solution, and obtain the entire branch we rerun
the continuation starting from the Kuramoto synchronization transition that Auto has detected:
The way Auto works, we cannot initialize the new simulation at `sym('LP1')`, we have to 
specify the original solution instead:

```
sol = run(new_branch('LP2'),DS="-")
```

Notice we had to change direction again with `DS="-"` (since we were going backwards).

Again, auto computes authomatically additional solutions from pitchfork bifurcations (we can avoid this turning off the detection of new branches).
From these two branches, we are only interested on the first, as we already have the second:

```
sym=sol[0]
plt.plot(asym['K'],asym['Ra'],'black')
plt.plot(asym['K'],asym['Rb'],'black')
plt.plot(sym['K'],sym['Ra'],'black')
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
```

We observe a new pitchfork (`BP`) of the homogeneous state at $K\approx 15.67$.
Auto did not check this branch. We can continue it with: 

```
another_branch=run(sol('BP2'),ISW=-1,STOP=['BP1'])
```

Notice the `ISW=-1` to force a branc switch, otherwise auto computes the solution we already know. Also we use `STOP=['BP2']` to avoid recomputing solutions already known for us.

This continuation provides two new branches. One is a branch of (unstable) assymetric fixed points that vanish on a subcritical pitchfork bifurcation to the antiphase state, which 
is a new solution for us: 

```
asym2 = another_branch[0]
antiphase = another_branch[1];
```

Let's visualize this:

```
plt.plot(asym['K'],asym['Ra'],'black')
plt.plot(asym['K'],asym['Rb'],'black')
plt.plot(sym['K'],sym['Ra'],'black')
plt.plot(asym2['K'],asym2['Ra'],'gray')
plt.plot(asym2['K'],asym2['Rb'],'gray')
plt.plot(antiphase['K'],antiphase['Ra'],'blue')
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
```

We should continue the antiphase solution downstream and this would have provided all the relevant fixed points in the system for this value of $p$. For instance, we can do

```
antiphase = run(another_branch('UZ1'),DS="-",ISW=1,STOP=[])
```

The plotting commands used previously should show now the full antiphase state.
Notice that here auto breaks when $R=0$ (as it should).
Also notice a Hopf bifurcation from the antiphase state. We save it for later use:

```
hopf2 = antiphase('HB1')
```

### Limit-cycles

Now, let's turn our attention back to the limit-cycle solutions
emerging from the Hopf of the chimera states. This Hopf bifurcation is at

```
hopf['K']
```

Auto can authomatically continue the resulting limit-cycles.
We just have to specify that we are interested in periodic orbits with `IPS=2`.
We also turn on the detection of bifurcations from periodic orbits with `ISP=2`:

```
lc=run(hopf,IPS=2,ISP=2,ICP=[1,11,3,4,5],NMX=50000,ISW=1,DSMAX=0.01,NTST=200,NCOL=7, STOP=['BP1'])
```

Now, the simulation halts at $K\approx 7.88$ without detecting any bifurcation.
However we can see that the period of the oscillation is quite large, indicating
a possible homoclinic bifurcation.

A closer inspection shows that the period orbit is colliding with the homogeneous state.
This can also be seen here if we plot the maxima and minima of the limit cycle,
together with the fixed points:

```
plt.plot(asym['K'],asym['Ra'],'black')
plt.plot(asym['K'],asym['Rb'],'black')
plt.plot(sym['K'],sym['Ra'],'black')
plt.plot(asym2['K'],asym2['Ra'],'gray')
plt.plot(asym2['K'],asym2['Rb'],'gray')
plt.plot(antiphase['K'],antiphase['Ra'],'blue')
plt.plot(lc['K'],lc['MAX Ra'],'r', lc['K'],lc['MIN Ra'],'r',lc['K'],lc['MAX Rb'],'r',lc['K'],lc['MIN Rb'],'r')
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
```

Let's go now to the Hopf bifurcation of the antiphase state, and see if it can complete the picture we have here:

```
lc2=run(hopf2,IPS=2,ISP=2,ICP=[1,11,3,4,5],NMX=50000,ISW=1,DSMAX=0.01,NTST=200,NCOL=7, STOP=['BP1'])
```

So, finally, we have that:

```
plt.plot(asym['K'],asym['Ra'],'black')
plt.plot(asym['K'],asym['Rb'],'black')
plt.plot(sym['K'],sym['Ra'],'black')
plt.plot(asym2['K'],asym2['Ra'],'gray')
plt.plot(asym2['K'],asym2['Rb'],'gray')
plt.plot(antiphase['K'],antiphase['Ra'],'blue')
plt.plot(lc['K'],lc['MAX Ra'],'r', lc['K'],lc['MIN Ra'],'r',lc['K'],lc['MAX Rb'],'r',lc['K'],lc['MIN Rb'],'r')
plt.plot(lc2['K'],lc2['MAX Ra'],'r', lc2['K'],lc2['MIN Ra'],'r',lc2['K'],lc2['MAX Rb'],'r',lc2['K'],lc2['MIN Rb'],'r')
plt.ylabel('Ra,Rb')
plt.xlabel('K')
plt.show()
```

We see the two limit-cycle branches are about to join, but the continuation from the antiphase solution breaks down when the orbit touches $R=0$. To obtain the trajectory in this narrow space one can perform simulations of the system. In spite of this (numerical) constrain, a detailed analysis shows that the limit-cycles from both Hopf solutions join at a Double Homoclinic bifurcation (see manuscript for more details).

All these solutions can be exported either using Python or directly with the auto commands, e.g., `save(lc,"lc")`. The computed data contains more information, such as the stability of each solution and the corresponding eigenvalues, see the official documentation for more information.

Auto generates several auxiliary files that might not be needed.
Before closing, do not forget to  `clean()`!


## Installing `auto-07p`

In order to install `auto-07p` in Linux from the official repository you can use:

```
mkdir auto-07p
git clone https://github.com/auto-07p/auto-07p auto-07p
cd auto-07p
./configure
make
make install
```

Depending on your system, there might be some conflicts.
I suggest installing a minimal version without the provided plotting tools,
as they require some dependencies that are outdated or conflict with current packages.
To do so, diable them in the configuration step of the previous instructions:

```
./configure --enable-plaut04=no --enable-plaut=no --enable-plaut-qt=no
```

Some other conflicts might appear anyhow, please see the official documentation.

If everything goes according to plan, typing `auto` in a terminal should start the interface.






