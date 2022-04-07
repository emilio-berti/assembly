
<!-- badges: start -->

[![R-CMD-check](https://github.com/emilio-berti/assembly/workflows/R-CMD-check/badge.svg)](https://github.com/emilio-berti/assembly/actions)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Example workflow

Following food web terminology, I talk here about *resources* and
*consumers*. In the case of plant-pollinator networks, this is analogous
as replacing *resource* with *plant* and *consumers* with *pollinator*.
In this case, however, some filtering steps may be unnecessary and
unwanted. **More about this in a separate section.**

I show an example workflow that is composed of the following steps:

1.  draw random species from a metaweb
2.  impose resource filtering, i.e. each basal species must be consumed
    and each consumer must feed on a resource.
3.  impose limiting similarity filtering, i.e. consumers are constantly
    replaced by others following a probability distribution proportional
    to their similarity of interactions.

## Draw random species from a metaweb

Loading the required libraries and set a random seed:

``` r
library(assembly)
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union

set.seed(1234)
```

Load the dataset *adirondack* that comes with *assembly*:

``` r
data(adirondack)
```

*adirondack* is the Adirondack Lakes metaweb as obtained from the
GATEWAy database.

Define the number of species for the local community:

``` r
S <- 50 #species richness
```

To draw a random community, use the function `draw_random_species()`:

``` r
sp <- draw_random_species(S, colnames(adirondack))
sum(colSums(adirondack[sp, sp]) == 0) #20 basal speciesshow_fw(sp, adirondack, title = "Random")
#> [1] 20
plot(graph_from_adjacency_matrix(adirondack[sp, sp]), vertex.label = NA)
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
#> [1] "Nephrocytium sp."      "Alona rectangula"      "Scenedesmus dimorphus"
#> [4] "Coelosphaerium sp."    "Scenedesmus serratus"  "Kirchneriella lunaris"
```

`assembly:::.find_replacements()` find suitable replacement for the
isolated species:

``` r
assembly:::.find_replacements(sp,
                              assembly:::.find_isolated(sp, adirondack),
                              adirondack,
                              keep.n.basal = TRUE)
#> [1] "Osmerus mordax"     "Cosmarium sp."      "Staurastrum sp."   
#> [4] "Desmidium sp."      "Anabaena sp."       "Xanthidium armatum"
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
#> [1] 6
```

Usually, a proper food web has only one component, i.e. all species are
connected by a path. Having more than one component means that the food
web is actually made of several disconnected communities. In the case
above, it also means that at least one of this disconnected communities
is composed of only one isolated species.

## Resouce filtering

To impose the resource filtering, I call the function
`resource_filtering()`. This takes as input the species names, the
metaweb, and an optional argument *keep.n.basal* to specify weather the
original number of basal species should be kept constant (default =
`FALSE`). **NOTE this may not be implemented correctly**

Behind the curtain, `resource_filtering()` calls the hidden functions as
a way to compress code and make it consistent. That’s why you shouldn’t
bother too much about hidden functions: they’re there because they’re
useful in the development of the package, rather than for your usage. If
they’re useful for you and you understand how they work, use them.

``` r
sp_resource <- resource_filtering(sp, adirondack, keep.n.basal = TRUE)
show_fw(sp_resource, adirondack, title = "Resource filtering")
```

<img src="man/figures/README-resource-1.png" width="50%" style="display: block; margin: auto;" />

``` r
plot(graph_from_adjacency_matrix(adirondack[sp_resource, sp_resource]),
     vertex.label = NA)
```

<img src="man/figures/README-resource-2.png" width="50%" style="display: block; margin: auto;" />

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
#>  [1] "Chydorus bicornutus"        "Diaptomus sicilis"         
#>  [3] "Lepomis gibbosus"           "Colletheca mutabilis"      
#>  [5] "Xanthidium sp."             "Diceras sp."               
#>  [7] "Trichotria tetractis"       "Conochiloides dossuarius"  
#>  [9] "Kelicottia longispina"      "Aphanothece sp."           
#> [11] "Ankistrodesmus sp."         "Daphnia pulex"             
#> [13] "Keratella cochlearis"       "Crucigenia rectangularis"  
#> [15] "Coregonus clupeaformis"     "Melosira sp."              
#> [17] "Cryptomonas ovata"          "Chrosomus eos"             
#> [19] "Eucyclops agilis"           "Ascomorpha ecaudis"        
#> [21] "Rhinichthys atratulus"      "Trichocerca pusilla"       
#> [23] "Lecane mira"                "Nitzschia sp."             
#> [25] "Sida crystallina"           "Rhizosolenia eriensis"     
#> [27] "Diaphanosoma birgei"        "Dictyosphaerium pulchellum"
#> [29] "nanoflagellates "           "Conochiloides unicornis"   
#> [31] "Polyarthra euryptera"       "Merismopedia tenuissima"   
#> [33] "Cyclops scutifer"           "Lecane sp."                
#> [35] "Notemigonus crysoleucas"    "fish eggs"                 
#> [37] "Keratella testudo"          "Lepadella cristata"        
#> [39] "Keratella taurocephala"     "Trichocerca cylindrica"    
#> [41] "fish fry"                   "Diaptomus leptopus"        
#> [43] "Scenedesmus arcuatus"       "Kelicottia bostoniensis"   
#> [45] "benthic detritus"           "Xanthidium armatum"        
#> [47] "Gomphonema sp."             "Peridinium wisconsinense"  
#> [49] "Arthrodesmus incus"         "Mesocyclops edax"
```

## Limiting similarity filtering

**To come**

## Example 1: Trophic levels in random and filtered communities

For 50 local communities I:

1.  draw random species
2.  apply the resource filtering procedure
3.  calculate the trophic level of the species in the local communities
4.  I calculate average and maximum trophic levels within each community

``` r
# TL_random <- matrix(NA, S, 50)
# TL_resource <- matrix(NA, S, 50)
# for (i in seq_len(50)) {
#   sp <- draw_random_species(S, colnames(adirondack))
#   sp_resource <- resource_filtering(sp, adirondack, keep.n.basal = TRUE)
#   TL_random[, i] <- ATNr::TroLev(adirondack[sp, sp])[, 1]
#   TL_resource[, i] <- ATNr::TroLev(adirondack[sp_resource, sp_resource])[, 1]
# }
# plot(colMeans(TL_random), colMeans(TL_resource),
#      xlim = c(1.1, 2.2), ylim = c(1.1, 2.2),
#      main = "Average trophic level",
#      xlab = "Random communities",
#      ylab = "Resource filtered communities")
# abline(0, 1)
# plot(apply(TL_random, 2, max), apply(TL_resource, 2, max),
#      xlim = c(2.2, 9.2), ylim = c(2.2, 9.2),
#      main = "Maximum trophic level",
#      xlab = "Random communities",
#      ylab = "Resource filtered communities")
# abline(0, 1)
```

Trophic level tend to be higher in resource-filtered communities
compared to random communities. This is due to isolated consumers
(trophic level = 1) in the random communities and that were replaced by
connected consumers when imposing resource filtering. It’s also evident
when computing the number of connected components:

``` r
show_fw(sp, adirondack)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="50%" style="display: block; margin: auto;" />

``` r
assembly:::.components(sp, adirondack)
#> [1] 6

show_fw(sp_resource, adirondack)
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="50%" style="display: block; margin: auto;" />

``` r
assembly:::.components(sp_resource, adirondack)
#> [1] 1
```
