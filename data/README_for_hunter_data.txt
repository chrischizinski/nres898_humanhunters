
Data set used in THE UNIQUE SPATIAL ECOLOGY OF HUMAN HUNTERS by
Atle Mysterud, Inger Maren Rivrud, Vetgard Gundersen, Christer M. Rolandsen and Hildegunn Viljugrein

Data provided by Atle Mysterud, email: <atle.mysterud@ibv.uio.no>.

See Material and Methods of the paper for detailed information about the data.

1. Spatial analyses run with R-script "SpatialAnalyses_The_unique_spatial_ecology_of_human_hunters.R"
1a. The file datS.txt is tab-separated and can be read into R by read.table(file="datS.txt",sep="\t",header=T)

Variables in columns are:
   Region: total number of disseminated Lyme borreliosis cases (includes Neuro_cases) reported in specific municipality and year
   fylkenr: running identity number of county
   idnr: running identity number of municipality (corresponding to the spesific row number of the adjacency matrix)
   densMoose16: moose density
   densRed16: red deer density
   densRoe16: roe deer density
   densreinhunters: reindeer density
   densmoosehunters: density of moose hunters
   densredhunters: density of red deer hunters
   densroehunters: density of roe deer hunters
   densHuman17: density of inhabitants
   kdensrein17: reindeer density
   reinhunters17: number of reindeer hunters
   moosehunters17: number of moose hunters
   redhunters17: number of red deer hunters
   roehunters: number of roe deer hunters
   humans2017: number of inhabitants
   totareal: municipality size

1b. The file adj.red.txt is tab-separated and can be read into R by read.table(file="adj.red.txt",sep="\t")
This is a 424X424 matrix defining whether each pair of municipalities is a neighbour (1), through shared border, or not (0).

2. Temporal analyses run with R-script "TemporalAnalyses_The_unique_spatial_ecology_of_hunters.R"
The file tdat.txt is tab-separated and can be read into R by read.table(file="tdat.txt",header=T,sep="\t")

Variables in columns are:
   Region: West, East, South and North
   County: county identity number for hunter home address
   year: specific year in period 2003-2017
   R_redhunters: change in number of red deer hunters for a given county from previous (t-1) to present (t) year
   R_moosehunters:change in number of moose hunters for a given county from previous (t-1) to present (t) year
   R_redharvest.West_tl1: change in red deer harvest in region West from previous to present year, lagged 1 year
   R_redharvest.Region_tl1: change in red deer harvest (sum region) from previous to present year, lagged 1 year


































