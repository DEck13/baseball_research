
# baseball_research

My research on baseball.  


### Challenging nostalgia

The Challenging_nostalgia directory contains a raw copy of my paper 
"Challenging nostalgia and performance metrics in baseball" which was recently 
published at CHANCE (see online-version.pdf).  This directory 
also includes a technical report that makes all of the calculations that were 
performed in the paper transparent.  This technical report is fully reproducible.  

To cite this paper:

  Eck, D.J. (2020). Challenging nostalgia and performance metrics in baseball.
  CHANCE, 33 (1), 16--25.

To cite the technical report:

  Eck, D.J. (2020).  Supporting data analysis for "Challenging nostalgia and 
  performance metrics in baseball."
  
A Shiny app that allows others to produce their own analyses within the framework of 
this paper can be viewed [here](https://deck13.shinyapps.io/challenging_baseball_nostalgia/).


### Era adjustment

[with Shen Yan] The era_adjustment directory contains tools which offer a way to properly adjust statistics of players from different eras under reasonable assumptions. The basic idea is to balance how well one did "vs their peers" against the size of their eligible MLB talent pool. The balancing of these two quantities provides more realistic "vs their peers" statistics, players who stood far above their peers get penalized if the era in which they did was sparsely populated. 

This directory is a work in progress.


### Spray distributions

[with Charlie Young and David Dalpiaz] The spraydists directory contains a  
manuscript on spray chart distributions. Spray chart distributions are 2-dimensional 
distributions of the final batted ball locations corresponding to a batter vs. pitcher matchup of 
interest. A shiny app that allows users to study batter vs. pitcher matchups visually and numerically can be viewed [here](https://seam.stat.illinois.edu). 

### Paradoxical relationship of fWAR and wins

The WinsvfWAR directory contains analyses which show that in any given season, a team's 
total fangraphs wins above replacement (fWAR) is a great single predictor for a team's 
total wins. We also show that a significant number of franchises consistently 
underperform or overperform their fWAR. These two facts are in direct constrast of each 
other and it is surprising that they can hold simultaneously. Data credit goes to 
[Fangraphs](https://www.fangraphs.com).


### Direct calculation of WAR

[with Sixian Li] The directWAR directory contains some exploratory analyses for a new 
simple approach to calculate a player's WAR that draws on ideas from causal inference. 
More to come.


### License

This research is under the GNU General Public License v3.0, this license does 
not include the file Gallupfavoritesport.png which appears in the 
Challenging_nostalgia directory.

