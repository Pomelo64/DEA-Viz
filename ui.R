library(shiny)
library(lpSolve)
library(Benchmarking)
library(smacof)
library(plotly)












# Define UI for DEA Viz application
shinyUI(
        navbarPage(
                title = 'DEA Visualization',
                
                
                tabPanel('Data upload', 
                         fluidPage(
                                 
                                 titlePanel("Upload Facotrs Data"),
                                 
                                 sidebarLayout(
                                         sidebarPanel(
                                                 helpText("Note:Please upload the inputs&outputs dataset in the",
                                                          "CSV format such that the inputs are the left-most columns",
                                                          "and the outputs are at the right-most columns of the dataset.",
                                                          "In the right side panel it is possible to check the format ",
                                                          "of the dataset."), 
                                                 tags$hr(),
                                                 
                                                 fileInput('factors_datafile', 'Upload the I&O Dataset',
                                                           accept=c('text/csv', 
                                                                    'text/comma-separated-values,text/plain', 
                                                                    '.csv')),
                                                 
                                                 numericInput("num_of_inputs", "Number of Input factors", 1), 
                                               
                                                 checkboxInput('header', 'Header', TRUE),
                                                 checkboxInput('dmu_labels', 'DMU Labels', FALSE),
                                                 
                                                 radioButtons('sep', 'Separator',
                                                              c(Comma=',',
                                                                Semicolon=';',
                                                                Tab='\t'),
                                                              ','),
                                                 
                                                 radioButtons('dec', 'Decimal Symbol',
                                                              c(Comma=',',
                                                                Dot='.'),
                                                              '.'),
                                                 
                                                 radioButtons('quote', 'Quote',
                                                              c(None='',
                                                                'Double Quote'='"',
                                                                'Single Quote'="'"),
                                                              '"'),
                                                 tags$hr(),
                                                 
                                                 
                                                 
                                                 helpText("Note: When the data shown in the right panel is what ",
                                                          "it is supposed to be, each variable in one column and",
                                                          "inputs are separated correctly from outputs,",
                                                          "then press the submit button."),
                                                          
                                                 actionButton("submit_button","Submit")
                                         ),
                                         
                                       
                                         mainPanel(
                                                 h4("Dataset Evaluation"),
                                                 verbatimTextOutput("dataset_evaluation_message"),
                                                 
                                                 h4("Dataset Description"),
                                                 verbatimTextOutput("factors_info"), 
                                                 
                                                 h4("Inputs Factors"),
                                                 tableOutput("inputs_table"),
                                                 
                                                 h4("Outputs Factors"),
                                                 tableOutput("outputs_table")
                                                 
                                                
                                         )
                                 )
                         )),
        
                navbarMenu('Plots',
                           tabPanel("CEM MDU",
                                    fluidPage(
                                            
                                            # Application title
                                            titlePanel("Cross-Efficiency Unfolding"),
                                            
                                            sidebarLayout(
                                                    
                                                    sidebarPanel(
                                                            #action button 
                                                            #size 
                                                            #transparency
                                                            #labels 
                                                            
                                                            helpText("Note: By pressing the 'plot' button, the uploaded data",
                                                                     "would be visualized on the right panel using cross-efficiency",
                                                                     "unfolding. Afterwards, the below widgets could help to have",
                                                                     "a better final map."), 
                                                            tags$hr(),
                                                            radioButtons("cem_approach",label = "CEM Approach", choices = list("Benevolent", "Aggressive"),selected = "Benevolent"),
                                                            actionButton("cem_mdu_button","Plot"), 
                                                            tags$hr(),
                                                            
                                                            sliderInput("cem_row_point_size", label = "Row Objects Point Size", min = 1 , max = 4 , value = 1),
                                                            sliderInput("cem_col_point_size", label = "Col Objects Point Size", min = 1 , max = 4 , value = 1), 
                                                            sliderInput("cem_row_transparency", label = "Row Objects Transparency", min = 0.1 , max = 1 , value = 0.5),
                                                            sliderInput("cem_col_transparency", label = "Col Objects Transparency", min = 0.1 , max = 1 , value = 0.5),
                                                            
                                                            checkboxInput('cem_unfolding_labels', 'Labels', TRUE)
                                                            
                                                            
                                                           
                                                    ),
                                                    
                                                    
                                                    # Show the caption, a summary of the dataset and an HTML 
                                                    # table with the requested number of observations
                                                    mainPanel(
                                                            
                                                            #plotlyOutput("cem_mdu_plot"),
                                                            plotlyOutput("cem_mdu_plot"),
                                                            verbatimTextOutput("cem_mdu_info"),
                                                            #tableOutput("cem_mdu_info")
                                                            downloadButton('download_cem_unfolding', 'Download the Plot')
                                                            
                                                            
                                                           
                                                    )
                                            )
                                    ) 
                                    ),
                           
                           
                           tabPanel("Co-Plot",
                                    fluidPage(
                                            
                                            # Application title
                                            titlePanel("Co-Plot"),
                                            
                                            sidebarLayout(
                                                    
                                                    sidebarPanel(
                                                            
                                                            helpText("Note: By pressing the 'plot' button, the uploaded data",
                                                                     "would be visualized on the right panel using co-plot ",
                                                                     "method.(Adler & Raveh, 2008) Further , the below",
                                                                     "widgets could help to have a better final map."), 
                                                            tags$hr(),
                                                            actionButton("coplot_button","Plot"), 
                                                            tags$hr(),
                                                            
                                                            sliderInput("coplot_point_size", label = "Point Size", min = 1 , max = 4 , value = 1),
                                                            sliderInput("coplot_vector_size", label = "Vector Size", min = 1 , max = 4 , value = 1),
                                                            sliderInput("coplot_point_transparency", label = "Point Transparency", min = 0.1 , max = 1 , value = 0.5),
                                                            sliderInput("coplot_vector_transparency", label = "Vector Transparency", min = 0.1 , max = 1 , value = 0.5),
                                                            sliderInput("coplot_vector_treshold", label = "Vector Treshold", min = 0.1 , max = 0.95 , value = 0.5),
                                                            
                                                            checkboxInput('coplot_labels', 'Labels', TRUE)
                                                            
                                                    ),
                                                    
                                                    
                                                    # Show the caption, a summary of the dataset and an HTML 
                                                    # table with the requested number of observations
                                                    mainPanel(
                                                            
                                                            plotOutput("coplot_plot"),
                                                            #plotlyOutput("cem_mdu_plot"),
                                                            verbatimTextOutput("coplot_info"),
                                                            #tableOutput("cem_mdu_info")
                                                            downloadButton('download_coplot', 'Download the Plot')
                                                            
                                                            
                                                            
                                                    )
                                            )
                                    ) 
                                    )
                         ),
                
                tabPanel('About')
        )
)