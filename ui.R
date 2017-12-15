
library(shiny)
library(DT)
library(shinycssloaders)
library(shinyjs)

ui <- tagList(
  useShinyjs(),
  navbarPage(
    #title  = actionLink("aboutLink", "LeafViz"),
    title = "sQTLviz",
    #title = a(id = "gitLink", href="https://github.com/davidaknowles/leafcutter/tree/master/leafviz","LeafViz", target = "_blank"), 
    id = "navBarPage",
    windowTitle = "sQTLviz",
    tabPanel("CommonMind Consortium sQTLs", 
             value = "resultsToPlot",
             # padding-top: 70px
             tags$style(type="text/css", "
#UCSC_button {
  margin-top: 10px;
  margin-bottom: 10px;
  display: block;
  text-align: center;
}
#UCSC {
  border-color: #ccc;
  color: #333;
}
             
#UCSC:hover {
  background-color: #E6E6E6;
}

#titlePanel {
  text-align: center;
}

#logo {
    margin: 0 auto;
    width: 4%;
    position: fixed;
    z-index: -1;
}

#title {
    size:huge;
}

")),
    tabPanel("Li et al TWAS SNPs", value = "yangResults"),
    tabPanel("Nalls et al PD GWAS SNPs", value = "GWASresults" ),
    tabPanel("About", value = "NA",
             
             fluidRow(id = "tabDiv",
                      column(
                        6,
                        offset=3,
                        div(
                          h2("What is this?"),
                          p( "This R", tags$a(href="https://shiny.rstudio.com","Shiny", target = "_blank"), 
                             "app presents and visualises splicing quantitative trait loci (sQTLs) discovered in 533 human frontal cortex samples.",
                             "Samples were provided by the Common Mind Consortium.",
                             "Data is available", a(href="https://www.synapse.org/#!Synapse:syn2759792", "here.", target = "_blank"),
                             "Reference: " 
                             ),
                          p(
                            "Fromer, M. et al. Gene expression elucidates functional impact of polygenic risk for schizophrenia.",em("Nat. Neurosci."), "10.1038/nn.4399 (2016)"
                          ),
                          p( "All 8944 associations between a SNP and a cluster of introns at FDR < 0.05 is presented.",
                             "For comparison, you can subset the table to retrieve:"
                          ),
                          tags$ul(
                             tags$li(strong("Li et al TWAS SNPs"), "SNPs from Figure 4 of Li et al"),
                             tags$li(strong("Nalls et al PD GWAS SNPs"), "SNPs that also appear in a recent meta-analysis of Parkinson's Disease genome-wide association studies (Nalls et al, 2014,", em("Nat Gen"), "doi: 10.1038/ng.3043)")
                          ),
                          h2("Methods"),
                          p( 
                            "Yang I Li, Garrett Wong, Jack Humphrey and Towfique Raj. Prioritizing Parkinson's Disease genes using population-scale transcriptomic data.", em("bioRxiv"), "https://doi.org/10.1101/231001 (2017)",
                            a(href="https://www.biorxiv.org/content/early/2017/12/08/231001", "Link", target="_blank" )
                          ),
                          h2( "Differential splicing events"),
                          p( "A cluster is defined as set of overlapping spliced junctions or introns.", 
                             "Clusters are initially ranked in the cluster results table by adjusted P value."
                          ),
                          tags$ul(
                            tags$li( strong("SNP -") , "the rs ID of the single nucleotide polymorphism (SNP)." ), 
                            tags$li( strong("Position -") , "the genomic coordinates (hg19) of the single nucleotide polymorphism (SNP)." ), 
                            tags$li( strong("Gene -"), "the HUGO gene sympbol for that gene."),
                            tags$li( strong("Cluster coordinates - "), "the coordinate span of the the cluster."  ),
                            tags$li( strong("q -"), "the Benjamini-Hochberg adjusted P value of the multinomial test of intron counts between conditions.")
                          ),
                          h2("Splicing event visualization"),
                          p("To view a cluster, click on a row in the cluster results table. This will start the plotting function.", 
                            "A cluster plot and table are generated each time a row in the cluster results is clicked."),
                          p("For a chosen cluster, the mean number of splice junctions supporting each intron is calculated for each genotype and then normalised as a fraction of the total counts.",
                            "Therefore for each condition the normalised counts will add up to 1.",
                            "Each intron is plotted as a line connecting its start and end coordinates with a thickness proportional to the displayed normalised count value.",
                            "The colour of the intron line indicates whether it is present in the annotation (red) or not (pink).",
                            "Any exons that are annotated as flanking or being contained within the cluster are added as rectangles to the plot.",
                            "If exons from multiple genes are connected by introns then their exons will be coloured according to their gene of origin."
                          ),
                          p("Each intron is presented as a row in the cluster view table."
                          ),
                          tags$ul(
                            tags$li( strong("chr, start, end"), "the genomic coordinates of the intron." ), 
                            tags$li( strong("verdict") , "the support given to two splice sites of the intron (start and end) by annotation." ),
                            tags$ul(
                              tags$li( em("annotated -"), "both splice sites are present in an annotated junction"),
                              tags$li( em("novel annotated pair -"), "both splice sites are annotated but are not annotated as being paired in a junction"),
                              tags$li( em("cryptic_fiveprime -"), "the 3\' splice site is annotated but the 5\' is not."),
                              tags$li( em("cryptic_threeprime -"), "the 5\' splice site is annotated but the 3\' is not ")
                            )
                          ),

                          h2("Gene-level visualization"),
                          p("This visualises all clusters discovered by Leafcutter that can be assigned to a particular gene.",
                            "Exons are taken from the provided annotation and plotted as black rectangles.",
                            "Each junction in each cluster is plotted as curved line with uniform thickness.",
                            "Note that the genomic coordinates are deliberately warped to give more space to the clusters."
                          ),
                          
                          h2("Junction-level visualization"),
                          p("The junction with the most significant association to the SNP is bolded in the cluster-level plot.",
                            "As well, the contribution of that junction to the cluster is calculated for each sample and split by genotype.",
                            "The Beta and q values come from the logistic regression."
                            ),
                          
                          h2("Acknowledgments"),
                          p("This work was supported by the US National Institutes of Health (NIH grant R01AG054005).", 
                            "We thank the patients and families who donated material for CommonMind Consortium data.", 
                            "The CommonMind Consortium data are available in the ",
                            a(href="https://www.synapse.org/#!Synapse:syn4923029", "CMC Knowledge Portal.", target = "_blank"),
                            "Data were generated as part of the CMC supported by funding from Takeda Pharmaceuticals Company Limited, F. Hoffman-La Roche Ltd and NIH grants R01MH085542, R01MH093725, P50MH066392, P50MH080405, R01MH097276, RO1-MH-075916, P50M096891, P50MH084053S1, R37MH057881 and R37MH057881S1, HHSN271201300031C, AG02219, AG05138 and MH06692.",
                            "Brain tissue for the study was obtained from the following brain bank collections: the Mount Sinai NIH Brain and Tissue Repository, the University of Pennsylvania Alzheimer's Disease Core Center, the University of Pittsburgh NeuroBioBank and Brain and Tissue Repositories and the NIMH Human Brain Collection Core."
                          ),
                          p("CMC Leadership:",
                            "Pamela Sklar, Joseph Buxbaum (Icahn School of Medicine at Mount Sinai), Bernie Devlin, David Lewis (University of Pittsburgh), Raquel Gur, Chang-Gyu Hahn (University of Pennsylvania), Keisuke Hirai, Hiroyoshi Toyoshiba (Takeda Pharmaceuticals Company Limited), Enrico Domenici, Laurent Essioux (F. Hoffman-La Roche Ltd), Lara Mangravite, Mette Peters (Sage Bionetworks), Thomas Lehner, Barbara Lipska (NIMH)."
                            )
                        )
                      )
             )
             
             ),
    
             
  
  # Application title
  div(id = "titlePanel",
    #HTML("<img id=logo src=squirtle.png>"),
    h1("sQTLviz: explore splicing QTLs", id = "titlePanel")
    ),
    
    
    fluidRow(
      column(8,offset=2,
        div(
          DT::dataTableOutput("all_clusters")
        )
      )
    ),
    withSpinner(div(
      column(7, offset = 1,
        h3(id = "title", "Whole gene visualization")
      ),
      column(12, offset = 0,
        plotOutput("select_gene_plot",width="100%", height = "350px") 
      ),
      fluidRow(
        column(7, offset = 1,
          div(
            h3(id = "title", "Cluster visualization"),
            h5(id = "subtitle", "Most significant junction is bolded")
          ),
          div(
           plotOutput("select_cluster_plot", width = "100%", height = "500px") 
          )
        ),
        column(3,
          h3(id = "title", "Most significant junction"),
          div(
            plotOutput("select_box_plot", width = "100%", height = "500px")
          )
        )
      ),
      fluidRow(
        column(8,offset=2,
          div(
        DT::dataTableOutput("junctionTable")
          )
        )
      ),
      div( id = "UCSC_button",
           htmlOutput("view_cluster_UCSC", inline = TRUE)
           )
    ), type = 8
    )
  )
)
