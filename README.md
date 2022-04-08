
<!-- badges: start -->

[![R-CMD-check](https://github.com/emilio-berti/assembly/workflows/R-CMD-check/badge.svg)](https://github.com/emilio-berti/assembly/actions)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Getting started

## Background

Following food web terminology, I talk here about *resources* and
*consumers*. In the case of plant-pollinator networks, this is analogous
as replacing *resource* with *plant* and *consumers* with *pollinator*.
In this case, however, some filtering steps may be unnecessary and
unwanted. **More about this in a separate section.** I also talk about
metawebs, which are the regional food web networks that emerge
considering all interactions that occur in local communities within the
region of interest.

The general scope of *assembly* is to simulate top-down assembly
processes. Top-down assembly means that all species are introduced into
a local community at once. The opposite, i.e. communities are built by
sequential introductions of one species, is called bottom-up assembly.
Bottom-up assembly is problematic for community composed of many
species, as the number of unique assembly sequences is
![S!](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%21 "S!"),
where
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S "S")
is the number of species in the metaweb. For 100 species, for instance,
there are
![\\sim 10^{157}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csim%2010%5E%7B157%7D "\sim 10^{157}")
unique sequences, making computations (and replications) virtually
impossible. Strikingly, Serván and Allesina (2021) showed that bottop-up
and top-down assembly are equivalent under some specific conditions.

Far from assuming this is the case in *assembly*, I simply want to
highlight that *assembly* only implements top-down assembly and that
none of the procedures in *assembly* can be considered steps in an
ecological sequence. Whenever I talk about *steps*, *moves*, and
*sequences* in *assembly*, I always refer to *procedures steps*,
*procedures moves*, and *procedures sequences*. Always keep that in mind
and don’t be fooled by the terminology; there is no bottom-up assembly
in *assembly*.

At this point, if you’re like me, you will ask *why not bottom-up
assembly?* The simple answer is: because it is computational unfeasible
for large communities. To make it feasible, it is tempting to restrict
the possible sequences to a (very) narrow space, only to 1,000,000
sequences
(![10^{-152}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;10%5E%7B-152%7D "10^{-152}")
of the total fraction of possible sequences for 100 species; much *much*
smaller than considering only a grain of sand simulating the whole
Earth). As this seems quite dangerous, potentially leading to biased
conclusions, I completed neglected bottom-up assembly processes.

## General workflow of *assembly*

In general, *assembly* implements this pipeline:

1.  Draw random species from a metaweb.
2.  Impose resource filtering, i.e. each basal species must be consumed
    and each consumer must feed on a resource.
3.  Impose limiting similarity filtering, i.e. consumers are replaced by
    others, depending on a probability distribution proportional to
    their similarity of interactions.

These steps are not necessarily sequential, i.e. limiting similarity can
be imposed independently from the resource filtering, but this is
programmatically harder. Because ecological communities are likely
filtered by resource availability before interaction competition, I
wrote *assembly* in a way that is straightforward to pass the output of
the resource filtering procedure to the limiting similarity one. It is,
however, possible to omit/invert procedures by performing extra-checks
on input/output data and implementing few pre-procedure steps. I will
show how to implement some of them at the end.

# Draw random species from a metaweb

Loading the required libraries and set a random seed

``` r
library(assembly)
#> Loading required package: igraph
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(igraph)

set.seed(123)
```

Load the dataset *adirondack* that comes with *assembly*:

``` r
data(adirondack)
```

*adirondack* is the Adirondack Lakes metaweb as obtained from the
GATEWAy database.

`draw_random_species(n, sp.names)` draws random species from a metaweb
and requires as input the number of species to draw (*n*) and the
species names (*sp.names*). If the metaweb does not have species names,
you can assign random ones with `name_metaweb(metaweb)`:

``` r
colnames(name_metaweb(matrix(as.numeric(runif(100) > .5), 10, 10)))
#>  [1] "ytvynywchp" "lyncngcwvz" "ouehsjrjlb" "jvltnqnvch" "nsoxqwkgow"
#>  [6] "zfngjefpxu" "wkdlnsygvz" "igbpmsxtog" "dahtypxvkp" "thcdtlvqjt"
```

If names are already present, this will throw an error:

``` r
tryCatch(name_metaweb(adirondack),
         error = function(e) print(e))
#> <simpleError in name_metaweb(adirondack): Names are already present>
```

When the metaweb has names, define the number of species for the local
community:

``` r
S <- 50 #species richness
```

And draw a random community, use the function `draw_random_species()`:

``` r
sp <- draw_random_species(S, colnames(adirondack))
show_graph(sp, adirondack)
```

<img src="man/figures/README-random-1.png" width="50%" style="display: block; margin: auto;" />

## Hidden functions

There are several hidden functions in *assembly*. The reason there are
hidden functions is because there is no need to call them directly.
Hidden functions can be accessed by prefixing the `assembly:::` (three
colon, not two). All hidden functions start with a dot `.`,
e.g. `assembly:::.basals()`.

In general, you should not be bothered by hidden functions and should
not call them directly, unless you have a good understanding of how they
operate. Nevertheless, I summarize them for clarity.

`assembly:::.basals()` get all basal species in the metaweb, and is
equivalent to subset the names of the metaweb where
`colSums(metaweb) == 0`:

``` r
identical(
  sort(intersect(assembly:::.basals(adirondack), sp)),
  sort(intersect(colnames(adirondack)[colSums(adirondack) == 0], sp))
)
#> [1] TRUE
```

`assembly:::.consumers()` and `assembly:::.top()` return the consumers
and top consumers of the metaweb, respectively.

`assembly:::.find_isolated()` returns the species that are isolated in
the local community:

``` r
assembly:::.find_isolated(sp, adirondack)
#> [1] "Halotheca sp."             "Conochiloides hippocrepis"
#> [3] "Keratella serrulata"       "Ceriodaphnia quadrangula" 
#> [5] "Polyarthra euryptera"      "Lepadella triptera"       
#> [7] "Prosopium cylindraceum"    "Catostomus catostomus"
```

`assembly:::.find_replacements()` find suitable replacement for the
isolated species:

``` r
assembly:::.find_replacements(sp,
                              assembly:::.find_isolated(sp, adirondack),
                              adirondack,
                              keep.n.basal = TRUE)
#> [1] "Ploesoma hudsoni"            "fish fry"                   
#> [3] "Salmo rutta"                 "Conochiloides unicornis"    
#> [5] "Keratella crassa"            "Coregonus artedii"          
#> [7] "Salvelinus fontinalis small" "Chroococcus sp."
```

If *keep.n.basal* is TRUE (default = FALSE), then the original number of
basal species will not change.

`assembly:::.move()` performs a move in the limiting similarity
procedure (more about this later):

``` r
tryCatch(assembly:::.move(sp, adirondack, t = 1),
         error = function(e) print(e))
#> <simpleError in assembly:::.move(sp, adirondack, t = 1): Isolated species detected in input>
```

This call to `assembly:::.move()` fails because isolated species are
detected in the input. This is a desired property of the function,
i.e. it fails when there is an unexpected behavior. All hidden functions
have some kind of behavior-check, which is a safety net to assure the
code is doing what you asked for.

Finally, `assembly:::.components()` returns the number of connected
components in the graph of the local community:

``` r
assembly:::.components(sp, adirondack)
#> [1] 5
```

Usually, a proper food web has only one component, i.e. all species are
connected by a path. Having more than one component means that the food
web is actually made of several disconnected communities. In the case
above, it also means that at least one of this disconnected communities
is composed of only one isolated species.

# Resouce filtering

To impose the resource filtering, call the function
`resource_filtering()`. This takes as input the species names
(*sp.names*), the metaweb (*metaweb*), and an optional argument
*keep.n.basal* to specify weather the original number of basal species
should be kept constant (default = `FALSE`). **NOTE this may not be
implemented correctly**

Behind the curtain, `resource_filtering()` calls the hidden functions as
a way to compress code and make it consistent. That’s why you shouldn’t
bother too much about hidden functions: they’re there because they’re
useful in the development of the package, rather than for your usage. If
they’re useful for you and you understand how they work, use them.

``` r
sp_resource <- resource_filtering(sp, adirondack, keep.n.basal = TRUE)
show_graph(sp_resource, adirondack)
```

<img src="man/figures/README-resource-1.png" width="50%" style="display: block; margin: auto;" />

Now the local community is fully connected, i.e. basal species always
have a consumer and consumers always have an available resource. It’s
possible to check this manually calling the hidden functions and working
on the adjacency matrix of the local community:

``` r
bas <- intersect(sp_resource, assembly:::.basals(adirondack))
cons <- intersect(sp_resource, assembly:::.consumers(adirondack))
all(rowSums(adirondack[bas, cons]) > 0)
#> [1] TRUE
all(colSums(adirondack[union(bas, cons), cons]) > 0)
#> [1] TRUE
```

Usually you don’t need to perform these checks, as I implemented them
within `resource_filtering()`. I also implemented a check for
disconnected components, to make sure that the resulting community has
no isolated species and only one actual community.

Bonus: because of these checks, now it is safe to perform a move of the
limiting similarity procedure:

``` r
assembly:::.move(sp_resource, adirondack, t = 1)
#>  [1] "Cosmarium sp."               "Alona costata"              
#>  [3] "Anabaena flos-aquae"         "Fragilaria sp."             
#>  [5] "Holopedium gibberum"         "Quadrigula closterioides"   
#>  [7] "Crucigenia quadrata"         "Chlamydomonas sp."          
#>  [9] "Chrysosphaerella longispina" "Euchlanis sp."              
#> [11] "Arthrodesmus octocornis"     "Euglena acus"               
#> [13] "Cyclops scutifer"            "Polyarthra major"           
#> [15] "Staurastrum megacanthum"     "Scenedesmus sp."            
#> [17] "Crucigenia crucifera"        "Bosmina longirostris"       
#> [19] "Ascomorpha ecaudis"          "Sida crystallina"           
#> [21] "Pediastrum tetras"           "Rhinichthys atratulus"      
#> [23] "Schroederia setigera"        "Dinobryon sp."              
#> [25] "Cocconeia sp."               "Euastrum sp."               
#> [27] "Fragilaria crotonensis"      "Sphaerocystis schroeteri"   
#> [29] "Peridinium wisconsinense"    "Oncorhynchus mykiss"        
#> [31] "Daphnia galeata"             "Micropterus salmoides"      
#> [33] "Peridinium limbatum"         "Microsystis sp."            
#> [35] "Pimephales promelas"         "Tropocyclops prasinus"      
#> [37] "Peridinium inconspicuum"     "Carteria sp."               
#> [39] "Tabellaria flocculosa"       "Aphanothece sp."            
#> [41] "Diatoma sp."                 "Micrasterias sp."           
#> [43] "Diaphanosoma birgei"         "Salmo rutta"                
#> [45] "Kelicottia bostoniensis"     "fish fry"                   
#> [47] "Ankistrodesmus spiralis"     "Diaptomus minutus"          
#> [49] "Colletheca mutabilis"        "Lepomis gibbosus"
```

# Limiting similarity filtering

The limiting similarity filtering is composed of a series of individuals
moves:

1.  A metric
    (![J](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;J "J"))
    representing the similarity of interaction is calculated for each
    species in the community.
2.  One species is removed with probability proportional to their
    ![J](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;J "J").
3.  The species removed is replaced by another species selected at
    random from the metaweb.

Each move can then be accepted of discarded based on a probabilistic
acceptance criterion, the Metropolis-Hasting algorithm (see below for
details). If the move is accepted, the new community will differ from
the starting one only in the removed/replaced species. If the move is
discarded, the new community is identical to the starting one. Each move
takes as input the community of the previous move (or the initial
community for the first move) and returns as output the new community
(or the original one).

The limiting similarity procedure is simply a series of moves (accepted
or not). The starting community for the limiting similarity procedure is
usually a community that has undergone already a resource filtering,
imposing thus that trophic competition comes after resource
availability. However, it is possible to skip the resource filtering
step; the only requirement of the starting community for the limiting
similarity filtering is that it does not have isolated species and has
only one connected components.

To impose the limiting similariy filtering, call the function
`similarity_filtering()`. This function requires the species names as
input (*sp.names*), the metaweb (*metaweb*), the argument *t* (default =
0), which is the temperature of the Metropolis-Hastings algorithm, and
*max.iter* (default = 1,000), which is the maximum number of moves
allowed.

``` r
sp_sim <- similarity_filtering(sp_resource, adirondack, t = 1, max.iter = 10)
show_graph(sp_sim, adirondack)
```

<img src="man/figures/README-limiting-1.png" width="50%" style="display: block; margin: auto;" />

## Metropolis-Hasting algorithm

To accept a move, it must pass a Metropolis-Hasting algorithm. If the
similarity of the new community is lower than the similarity of the old
one, the move is always accepted. When the new similarity is higher than
the old similarity, the move is accepted if:

![
e^{\\left( 1 - \\frac{similarity\_{new}}{similarity\_{old}} \\right) \\frac{1}{t}} &gt; \\mathcal{U}(0, 1)
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ae%5E%7B%5Cleft%28%201%20-%20%5Cfrac%7Bsimilarity_%7Bnew%7D%7D%7Bsimilarity_%7Bold%7D%7D%20%5Cright%29%20%5Cfrac%7B1%7D%7Bt%7D%7D%20%3E%20%5Cmathcal%7BU%7D%280%2C%201%29%0A "
e^{\left( 1 - \frac{similarity_{new}}{similarity_{old}} \right) \frac{1}{t}} > \mathcal{U}(0, 1)
")

This means that, even when the new similarity is higher than the old one
and the new community has species with increased similarity of
interaction, it can still be accepted as a valid move based on a
probability density function. The probability of acceptance depends on
how much the new similarity is higher than the old one and by the
temperature parameter *t*. For increasing *t*, it is more likely to
accept a non-favorable move:

``` r
temp <- 10 ^ seq(-2, 1, by = .1)
ratio <- 10 ^ seq(-1, 2, by = .1)

move <- rep(NA, length(temp) * length(ratio))
i <- 1
for (t in temp){
  for (x in ratio) {
    move[i] <- metropolis.hastings(1, x, t)
    i <- i + 1
  }
}

d <- data.frame(Temp = rep(temp, each = length(ratio)),
                Ratio = rep(ratio, length(temp)),
                Move = as.numeric(move))
cols <- colorRampPalette(c("steelblue", "tomato"))
cols <- cols(length(unique(d$Temp)))
pal <- adjustcolor(cols, alpha.f = .2)
pal <- pal[sapply(d$Temp, \(x) which(unique(d$Temp) == x))]
plot(log10(d$Ratio), jitter(d$Move, factor = .2), 
     xlab = "Similiarity ratio (new / old)",
     ylab = "Accepted move",
     col = pal, pch = 20, frame = FALSE)
for (t in unique(d$Temp)) {
  x <- log10(d$Ratio[d$Temp == t])
  y <- d$Move[d$Temp == t]
  fit <- loess(y ~ x)
  lines(fit$x, fit$fitted, col = cols[which(unique(d$Temp) == t)])
}
legend(-1, .7, legend = seq(-2, 1, by = .5), 
       fill = colorRampPalette(c("steelblue", "tomato"))(7),
       title = "log10(t)")
```

<img src="man/figures/README-metro-1.png" width="80%" style="display: block; margin: auto;" />

The reason to include this algorithm is to avoid a purely deterministic
procedure and include some stochasticity in the process. However, if
this is unwanted, it can be removed (and the process made purely
deterministic), by specifying a very high temperature parameter *t*:

``` r
table(sapply(seq_len(1000), \(x) metropolis.hastings(1, 1e3, t = 1e9)))
#> 
#> TRUE 
#> 1000
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-servan2021tractable" class="csl-entry">

Serván, Carlos A, and Stefano Allesina. 2021. “Tractable Models of
Ecological Assembly.” *Ecology Letters* 24 (5): 1029–37.

</div>

</div>
