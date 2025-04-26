---
---

## Notes

- The Masuda paper supplementary actually includes an implementation of nMGA (which they ran to compare with their Laplace method)

- A reminder that the Boguna paper uses a linear approximation in order to compute the rates forward in time by delta t.
- that the Laplace GIllepsie algorithm is a double-draw algorithm: once as usual to determine the time as a function of rate, the other to draw the (Exponential) rate as a function of the actual distribution.
- That Grossman forgoes this by using rejection sampling to correct the distribution, instead of exactly comupting the hazard at t+delta t.

We show that the rejection method can be combined with exact delay methods: we will need to also show that there are competing comuptational tradeoffs for each process, which can be avoided by using each of the different two methods we use.




## Implementation Notes

- Masuda code has both NGMA and LGA code, we can use that to simplify our process. Problem is that they are hard-coded C(++?) code which is going to slightly diffcult to dientangle.
- We need to create a base Simulator class that sets up the basic initialisation and run loop signatures, and let each simulator implement its own .sim() methods (or abstract templates)

- This will allow us to more easily plug-and-play different simulators for comparison.

- A running concern is this delta-function hazard for delayed events. I dont think these "smooth" (L) and "linear" (nMGA) do a great approximating them, though we can try with a very narrow uniform distribution instead. (I know that rejection should be able to handle uniform, but its likely that nothing convolves to a uniform (and linear hazard approx to uniform kinda only works if dt is tiny))

- LGA
    - problem here is that since we draw the process rate over and over agai, the dict consolidation methods does not work as well
        - we have less repeats,
        - and we have a problem where we cannot consolidate multiple events, since we probably want to update their rates independently to save computation time.
- NMGA
    - need to work out how OFTEN these rates update
        - we run into the same problem that since the rates are varaible, we cannot consolidate with a single node, since we will not know how to split them.

-> Looks like the obivous way to construct a rate arr-ish thing with a tuple of (node, event_class) as the mapping key. We will also have to be aggressive with how we choose to delete keys from the dict -> will this be effective in Python?
    - with this, does this mean that the ratedict is overkill?
    - we can just walk a rate arr, extract the index, and then look up the key?

- we can construct a different scaffold for LGA; and for nMGA: maybe we need to also construct a scaffold for the functional aspect?? => 

i think nMGA works by computing a lnear approximation about delta T?
so we make the approx that

Surv(t) ~ Surv(t0) + haz(t0) (t - t0)

And then we need to solve Sum_i Surv_i(t) ~ [Uniform(0,1)]
This is pretty simple linear solve, but requires evaluating the Surv and haz fns at each t0 (i.e) at each 
we have a gamma hazard, so does this become difficult to evaluate? We also have to do this per reaction...? Does this make sense in the large N limit? This seems prohibitively expensive...
WE can mix this with the constant hazard processes actually.
Like have a static rate dict that we draw from like in GMax, and also a dynamic one that needs to be evaluated per event.
We have E -> I as our gamma dist, i cant imagine how poorly behaved this would be when we have a large N state (like if S -> T on a random uniform dist...)
Also nMGA only really works when dT is small, i.e. things are fast, and when we have to re-eval over and over again.

NMGA:
we linearly approximate the surivival function. The total rate becomes 

  exp [ - Sum_j log [ Surv_j(t_j) / Surv_j(t_j + dt) ]]
= exp [ - SUm_j log [ Surv_j(t_j) / [Surv_j(t_j) + haz_j(t_j) dt ]]]