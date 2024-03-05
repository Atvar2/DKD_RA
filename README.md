# Single-cell insights into the alleviation of diabetic nephropathy by Rosmarinic acid
Collection of R codes and horse perl scripts to perfrom analysis and data visualization in research article "Single-cell transcriptomics reveals the ameliorative effect of Rosmarinic acid on diabetic nephropathy-induced kidney injury by modulating oxidative stress and inflammation".

|         |                                                                  |
| ------- | ---------------------------------------------------------------- |
| Authors | Junhui Chen ([chenjunhui](https://github.com/Atvar2))         |
| Email   | <chenjhbio@163.com>                                           |
| License | [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)               |

## Citations
If you used the codes or reference part of the code in your research, please kindly cited following paper:
Chen et al. ABSP (2022) Single-cell transcriptomics reveals the ameliorative effect of Rosmarinic acid on diabetic nephropathy-induced kidney injury by modulating oxidative stress and inflammation

## Dependece packages
Following are a list of R packages that are used by the analysis pipeline. Before implement the codes, please comfirmed the packages have been installed on your R platform.

- Seaurat version 4.1.0
- DoubletFinder version 2.02
- clusterProfiler version 4.0.0
- SCP version: 0.5.6 
- ggplot2 version: 3.3.5

The other R packages such as dplyr, should depend on R environment variable space to install the responding packages as following commands:
```bash
install.packages(dplyr) 
or
BiocManager::install(dplyr)
```
## Installation
Just download the codes and implement in R platform under linux or window system. For packages installation, mainly intallated through [conda](https://docs.conda.io/en/latest/) or build-in funciont of R install.packages() 
