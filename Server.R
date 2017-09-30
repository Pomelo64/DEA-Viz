library(shiny)
library(readr)
library(lpSolve)
library(Benchmarking)
library(smacof)
library(ggplot2)
library(ggrepel)
library(plotly)

set.seed(7)
#####
##CEM MDU Functions

all_OFs = function(input_mat, output_mat) {
        #this function receives the input and output levels and generates a matrix in which each 
        #row is the objective function coefficients of the corresponding unit. 
        #in other words, all_OFs[i,] is the OF coefficients of unit i
        #ready to be used in lp()
        
        # I assume that the row size of the input_mat and output_mat are equal
        # Although, it is possible to put a check here 
        input_mat = as.matrix(input_mat)
        output_mat = as.matrix(output_mat)
        
        number_of_units = nrow(input_mat)
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        
        #creating the final matrix 
        all_OFs = matrix(nrow = number_of_units, ncol =  number_of_outputs+number_of_inputs  )
        
        #filling the final matrix
        for (unit in 1:number_of_units) {
                all_OFs[unit,] = c(output_mat[unit,],rep(0,number_of_inputs))
                
        }
        
        all_OFs
        #now by picking every row
        #it is possible to add input and output names to this matrix. Maybe in next versions
        
}

A_mat = function(unit, input_mat, output_mat) {
        #this function returns A matrix - matrix of constraints - for DMUunit formulation. 
        # unit can be a number from 1 to nrow(input_mat)
        # This function can be developed in a way that it returns all A matrices for
        # all units in a data structure such as list. It would be faster but 
        # pre-mature optimization is the root of all devils. 
        
        input_mat = as.matrix(input_mat)
        output_mat = as.matrix(output_mat)
        
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        number_of_units = nrow(input_mat)
        #assuming that the nrow() of input_mat and output_mat are equal 
        
        A = matrix(nrow = number_of_units, ncol= number_of_inputs + number_of_outputs)
        
        #first constraint is always the numerator constraint
        A[1,]= c(rep(0,number_of_outputs),input_mat[unit,])
        #print(dim(A))
        input_mat = input_mat[-unit, ]
        output_mat = output_mat[-unit, ]
        
        d = number_of_units-1
        
        #print(d)
        for (i in 1:d) { 
                #print(i)
                A[(i+1), ]=c(output_mat[i, ],-1*input_mat[i, ])
        }
        
        A
        
}

constraints_directions = function(input_mat) {
        #the = and <= directions of A matrix 
        # assuming that nrow(input_mat)==nrow(output_mat)
        number_of_units = nrow(input_mat)
        
        t = c("==",rep("<=",number_of_units-1))
        
        t
        
}

rhs = function(input_mat){
        #returns the rhs values
        #assuming nrow(input_mat)==ncol(input_mat)
        
        number_of_units = nrow(input_mat)
        
        t=c(1,rep(0,number_of_units-1))
        
        t
}

simple_efficiency = function(input_mat,output_mat,epsilon_value = 0.00001){
        #this function is supposed to return the simple efficiency scores of all units. 
        #assuming that nrow(input_mat)==nrow(output_mat)
        number_of_units = nrow(input_mat)
        OF_coefficients = all_OFs(input_mat, output_mat)
        
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        d = number_of_inputs + number_of_outputs
        
        #A_lower = epsilons_mat(input_mat=input_mat,output_mat=output_mat,epsilon_value = 0.00001)
        A_lower = diag(d)
        
        C_upper =constraints_directions(input_mat = input_mat )
        C = c(C_upper,rep(">=",d))
        
        rhs_upper=rhs(input_mat = input_mat)
        rhs = c(rhs_upper, rep(epsilon_value,d))
        
        simple_eff = vector(length = number_of_units)
        
        for (i in 1:number_of_units){
                
                A_upper=A_mat(unit=i, input_mat=input_mat, output_map=output_mat)
                A = rbind(A_upper,A_lower)
                
                
                
                
                DEA_LP = lp(direction = "max",objective.in = OF_coefficients[i,],const.mat = A, const.dir = C, const.rhs = rhs  )
                #model$solution & $objval
                
                simple_eff[i] = DEA_LP$objval 
                
        }
        
        simple_eff  
        
        
}

CEM_unit = function(input_mat , output_mat , unit , epsilon = 0.00001){
        #this function must return the the benevolent optimum weights for the given unit 
        
        #requires Benchmarking library
        #requires lpSolve library
        
        #require(Benchmarking)
        require(lpSolve)
        
        number_of_units = nrow(input_mat)
        input_mat = as.matrix(input_mat)
        output_mat = as.matrix(output_mat)
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        
        #we need the simple efficiency of the given unit 
        #require(Benchmarking)
        #library(Benchmarking)
        
        # 7 October - In order to test the new function dea_4cem instead of dea() 
        #simple_eff = dea(X = input_mat , Y = output_mat , RTS = "crs", ORIENTATION = "in")
        #unit_simple_eff = simple_eff$eff[unit]
        eff_weight_mat = dea_4cem(input_mat = input_mat, output_mat = output_mat )
        simple_eff = eff_weight_mat[,1]
        unit_simple_eff = simple_eff[unit]
        
        
        output_mat_refined = output_mat[-unit,]
        OF = apply(output_mat_refined,2,sum)
        OF = c(OF,rep(0,number_of_inputs))
        
        d = number_of_inputs + number_of_outputs
        
        #preparing the A matrix
        A_middle = diag(d)
        A_upper=A_mat(unit=unit, input_mat=input_mat, output_mat=output_mat)
        A_last = c(output_mat[unit,],-1*unit_simple_eff*input_mat[unit,])
        A = rbind(A_upper,A_middle,A_last)
        #two new lines
        input_mat_refined = input_mat[-unit,]
        A[1,] = c(rep(0,number_of_outputs),apply(input_mat_refined,2,sum))
        
        #preparing the constraint directions
        C_upper = c("==",rep("<=",number_of_units-1))
        C_middle = rep(">=",d)
        C_last = "=="
        C = c(C_upper,C_middle,C_last)
        
        #preparing the RHS values
        rhs_upper = c(1,rep(0,number_of_units-1))
        rhs_middle = rep(epsilon,d)
        rhs_last = 0
        rhs_total = c(rhs_upper, rhs_middle,rhs_last)
        ## changed rhs name to rhs_total
        #print(cbind(A,C,rhs))
        
        #changed rhs name to rhs_total
        t = lp(direction = "max", objective.in = OF , const.mat = A , const.dir = C , const.rhs = rhs_total)
        #CEM_weights_of_unit = t$solution
        #just made the above line inactive, why not returning the solutions directly? 
        t$solution
}

CEM = function(input_mat, output_mat) {
        #this function returns the benevolent CEM 
        
        require(lpSolve)
        number_of_units = nrow(input_mat)
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        d = number_of_inputs+number_of_outputs
        
        CEM_opt_weights = matrix(nrow = number_of_units, ncol = d )
        
        for (unit in 1:number_of_units) {
                CEM_opt_weights[unit,] = CEM_unit(input_mat , output_mat , unit , epsilon = 0)
        }
        #print(CEM_opt_weights)
        CEM = matrix (nrow = number_of_units, ncol = number_of_units)
        
        #In the 35 Chinesse Cities dataset, CEM_unit() could not find feasible solution for 
        #the unit 4. In other words, the benevolent formulation of the unit4 was not feasible!
        #why? I don't know now (3-oct-2016). As the result, the dea() function returns infeasible 
        #solution with ZERO as all weights. So here if I detect ALL ZERO, I replace it with the optimum weights 
        #that dea() function returns for the problematic unit. These weights - which may not be very benevolent
        #will be used in generation of cross-efficiency matrix. 
        
        eff_weight_mat = dea_4cem(input_mat, output_mat) 
        for (row in 1:number_of_units){
                if ( sum(CEM_opt_weights[row,])==0 ) { 
                        
                        #temporary =dea(X = input_mat, Y = output_mat , RTS = "crs" , DUAL = TRUE, ORIENTATION = "in")
                        # Orientation has been added on 7 october, trying to resolve the bug of colombian hospitals
                        #CEM_opt_weights[row,] = c(temporary$vy[row,] , temporary$ux[row,])
                        CEM_opt_weights[row,] = eff_weight_mat[row,-1]
                }
        }
        
        
        for (row in 1:number_of_units){
                w_outputs = CEM_opt_weights[row,1:number_of_outputs]
                w_inputs = CEM_opt_weights[row,(number_of_outputs+1):d]
                #print(U)
                #print(V)
                #print("----")
                
                ## This can be heavily optimized with matrix operations rather than scalar
                CEM[row,] = apply(X = t(t(output_mat) * w_outputs), MARGIN = 1 , FUN = sum ) / apply(X = t(t(input_mat) * w_inputs), MARGIN = 1 , FUN = sum )
                
                #for (col in 1:number_of_units){
                #        CEM[row,col] = sum(w_outputs*output_mat[col,])/sum(w_inputs*input_mat[col,])
                #        #CEM[row,col] = sum(round(U*output_mat[col,],4))/sum(round(V*input_mat[col,],4))
                #}
                
        }
        
        CEM
        
}

dea_4cem = function(input_mat , output_mat) {
        # this function is supposed to replace dea() of Benchmarking
        # This function supposed to give back the optimum [ or a optimum set of] weights
        
        # since dea() is giving bulshit results for Colombian hospital case
        # indeed, the summation of the VIs is not equal to one! 
        # so it means the formulation that dea() is using is not the one that I need
        
        
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        #assuming that input_mat and output_mat have the same number of rows!
        number_of_units = nrow(input_mat)
        d = number_of_inputs + number_of_outputs
        
        eff_weight_mat = matrix(nrow = number_of_units, ncol = (d+1))
        
        #A_epsilon=(diag(d))
        
        for (unit in 1:number_of_units){
                
                OF = unlist(c(as.vector(output_mat[unit,]),rep(0,number_of_inputs)))
                
                A_upper = unlist(c(rep(0,number_of_outputs),input_mat[unit,]))
                A_middle = cbind(output_mat[-unit,],-input_mat[-unit,])
                A_bottom = OF 
                A_upper = unname(A_upper)
                A_middle = unname(A_middle)
                A_bottom = unname(A_bottom)
                #A_epsilon= unname(A_epsilon)
                #colnames(A_epsilon) = colnames(A_upper)
                A = as.matrix(rbind(A_upper, A_middle, A_bottom))
                #A = rbind(A, A_epsilon)
                
                #const_directions = c("==",rep("<=",(number_of_units)),rep(">=",d))
                const_directions = c("==",rep("<=",(number_of_units)))
                # epsilon = 0.0000001
                #RHS_values = c(1,rep(0,(number_of_units-1)),1,rep(0.0000001,d))
                RHS_values = c(1,rep(0,(number_of_units-1)),1)
                
                lp_model=lp(direction = "max", objective.in = OF , const.mat = A , const.dir = const_directions, const.rhs = RHS_values)
                
                eff_weight_mat[unit,] = c(lp_model$objval,lp_model$solution)
                
                
                
        }
        
        eff_weight_mat
        
        
        
}

CEM_unit_agg =  function(input_mat , output_mat , unit , epsilon = 0.00001){
        #this function must return the the aggresive optimum weights for the given unit 
        
        #requires Benchmarking library
        #requires lpSolve library
        
        require(Benchmarking)
        require(lpSolve)
        
        number_of_units = nrow(input_mat)
        input_mat = as.matrix(input_mat)
        output_mat = as.matrix(output_mat)
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        
        #we need the simple efficiency of the given unit 
        #require(Benchmarking)
        #library(Benchmarking)
        #simple_eff = dea(X = input_mat , Y = output_mat , RTS = "crs", ORIENTATION = "in")
        #unit_simple_eff = simple_eff$eff[unit]
        eff_weight_mat = dea_4cem(input_mat = input_mat, output_mat = output_mat )
        simple_eff = eff_weight_mat[,1]
        unit_simple_eff = simple_eff[unit]
        
        output_mat_refined = output_mat[-unit,]
        OF = apply(output_mat_refined,2,sum)
        OF = c(OF,rep(0,number_of_inputs))
        
        d = number_of_inputs + number_of_outputs
        
        #preparing the A matrix
        A_middle = diag(d)
        A_upper=A_mat(unit=unit, input_mat=input_mat, output_mat=output_mat)
        A_last = c(output_mat[unit,],-1*unit_simple_eff*input_mat[unit,])
        A = rbind(A_upper,A_middle,A_last)
        #two new lines
        input_mat_refined = input_mat[-unit,]
        A[1,] = c(rep(0,number_of_outputs),apply(input_mat_refined,2,sum))
        
        #preparing the constraint directions
        C_upper = c("==",rep("<=",number_of_units-1))
        C_middle = rep(">=",d)
        C_last = "=="
        C = c(C_upper,C_middle,C_last)
        
        #preparing the RHS values
        rhs_upper = c(1,rep(0,number_of_units-1))
        rhs_middle = rep(epsilon,d)
        rhs_last = 0
        rhs_total = c(rhs_upper, rhs_middle,rhs_last)
        ## changed rhs name to rhs_total
        #print(cbind(A,C,rhs))
        
        #changed rhs name to rhs_total
        t = lp(direction = "min", objective.in = OF , const.mat = A , const.dir = C , const.rhs = rhs_total)
        #CEM_weights_of_unit = t$solution
        #just made the above line inactive, why not returning the solutions directly? 
        t$solution
}

CEM_agg = function(input_mat, output_mat) {
        #this function returns the benevolent CEM 
        
        number_of_units = nrow(input_mat)
        number_of_inputs = ncol(input_mat)
        number_of_outputs = ncol(output_mat)
        d = number_of_inputs+number_of_outputs
        
        CEM_opt_weights = matrix(nrow = number_of_units, ncol = d )
        
        for (unit in 1:number_of_units) {
                CEM_opt_weights[unit,] = CEM_unit_agg(input_mat , output_mat , unit , epsilon = 0)
        }
        #print(CEM_opt_weights)
        CEM = matrix (nrow = number_of_units, ncol = number_of_units)
        
      
        
        for (row in 1:number_of_units){
                if ( sum(CEM_opt_weights[row,])==0 ) { 
                        temporary =dea(X = input_mat, Y = output_mat , RTS = "crs" , DUAL = TRUE)
                        
                        CEM_opt_weights[row,] = c(temporary$vy[row,] , temporary$ux[row,])
                }
        }
        
        
        for (row in 1:number_of_units){
                U = CEM_opt_weights[row,1:number_of_outputs]
                V = CEM_opt_weights[row,(number_of_outputs+1):d]
                #print(U)
                #print(V)
                #print("----")
                
                for (col in 1:number_of_units){
                        CEM[row,col] = sum(U*output_mat[col,])/sum(V*input_mat[col,])
                }
                
        }
        
        CEM
        
}


#####


shinyServer(function(input, output) {
        
        ##### 
        ##First reactive expression of the uploaded file 
        #to make the uploaded file a reactive expression
        datafile <- reactive({
                inFile <- input$factors_datafile
                
                if (is.null(inFile))
                        return(NULL)
                
               read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                                quote=input$quote, stringsAsFactors = FALSE, dec = input$dec)
                
                #read_delim(file = inFile$datapath, col_names = input$header , quote = input$quote,delim = input$sep  )
        })
        
        #####
        ## for the upload panel 
        
        # for generating the description of the dataset
        output$factors_info <- renderText({
                t<- datafile()
                
                paste("The dataset of",nrow(t),"DMUs, composed of",as.numeric(input$num_of_inputs),"inputs, and",ncol(t)-as.numeric(input$num_of_inputs),"outputs.")
        })
        
        # for rendering reactive table of inputs
        output$inputs_table <-renderTable({
                
                # inFile <- input$factors_datafile
                #if (is.null(inFile))
                #        return(NULL)
                #read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                #         quote=input$quote)[1:6,1:as.numeric(input$num_of_inputs)]
  
                datafile()[1:6,1:as.numeric(input$num_of_inputs)]
                
                
        })
     
        #for rendering reactive table of outputs    
        output$outputs_table <- renderTable({
               #inFile <- input$factors_datafile
                #if (is.null(inFile))
                #        return(NULL)
                #t <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                #         quote=input$quote)[1:6,]

                
                datafile()[1:6,(as.numeric(input$num_of_inputs)+1):ncol(datafile())]
        })
        
        ##### 
        ## Evaluation of the uploaded data and passing data to computation part 
        
        #returns two messages based on whether the dataset meets all the requirements or not
        dataset_evaluation <- eventReactive(input$submit_button, {
                
                
                #t <- datafile()
                t <- data.frame(sapply(datafile(),parse_number))
                dataset_evaluation_result <- vector(length = 3)
                names(dataset_evaluation_result) <- c("Numerical Factors","Non-negative Factors","No Missing Value")
                
                
                dataset_evaluation_result[1] <- ( sum(sapply(t,is.numeric)) == ncol(t) ) 
                dataset_evaluation_result[2] <- ( !(FALSE %in% (t>=0)) ) 
                dataset_evaluation_result[3] <- ( !(TRUE %in% is.na(t)) )
                
                ifelse(test = (FALSE %in% dataset_evaluation_result) , 
                       yes = "Error! There is something wrong with the dataset. It maybe either having non-numeric data, negative values, or missing values." ,
                       no = "Great! The dataset meets the requirements." )  
                
        })
        
        #sends the evaluation result as a message to UI
        output$dataset_evaluation_message <- renderText({
                
                dataset_evaluation()
                #t <- datafile()
                #t <- data.frame(sapply(t,parse_number))
                #sapply(t,class)
        })
        
        # send the verified dataset into a reactive expressing, so it can be used later on 
        final_dataset_reactive <- reactive({
                data.frame(sapply(datafile(),parse_number))
                #ifelse(test = dataset_evaluation() == "Great! The dataset meets the requirements.",yes = data.frame(sapply(datafile(),parse_number)) , no = NULL ) 
        })

        
# CEM MDU        
#####


        #generating benevolent CEM by pressing the 'cem_mdu_button' 
        cem_reactive <- eventReactive(input$cem_mdu_button, {
                
                number_of_inputs <- as.numeric(input$num_of_inputs)
                t<- final_dataset_reactive()
                
                switch(EXPR = input$cem_approach ,
                      "Benevolent" = CEM(input_mat = as.matrix(t[,1:number_of_inputs]) ,output_mat = as.matrix(t[,(number_of_inputs+1):ncol(t)] )) , 
                      "Aggressive" = CEM_agg(input_mat = as.matrix(t[,1:number_of_inputs]) ,output_mat = as.matrix(t[,(number_of_inputs+1):ncol(t)] )) )
                #CEM(input_mat = as.matrix(t[,1:number_of_inputs]) ,output_mat = as.matrix(t[,(number_of_inputs+1):ncol(t)] )) 
                
                 
        })
        
        
        
        #CEM MDU plot
        cem_unfolding <- reactive({
                t <- cem_reactive()
                
                 unfolding(delta = round((1-t),2),ndim = 2)
        })
        
        cem_unfolding_plot <- reactive({
                
                t <- cem_unfolding()
                g <- ggplot() + geom_point(data = data.frame(t$conf.row), aes(x = D1 , y = D2), color = "blue", size = input$cem_row_point_size, shape = 19, alpha = input$cem_row_transparency) + 
                        geom_point(data = data.frame(t$conf.col), aes(x = D1 , y = D2), color = "red",size = input$cem_col_point_size, shape = 24, alpha = input$cem_col_transparency) +
                        #geom_text_repel(data = data.frame(t$conf.row), aes(x = D1 , y = D2), color = "blue", alpha = 0.5 , label = label)+
                        #geom_text_repel(data = data.frame(t$conf.col), aes(x = D1 , y = D2), color = "orangered", alpha = 0.5 , label = label) +
                        coord_fixed(ratio = 1,  expand = TRUE)
                ggplotly(g)
                
        })
        
        output$cem_mdu_plot <- renderPlotly({
                
               # t <- cem_reactive()
                #label = 1:nrow(t)
                #t<- unfolding(delta = round((1-t),2),ndim = 2)
                
                cem_unfolding_plot()
                #ggplotly(g)
               
        })
        
        #Info of CEM unfolding
        output$cem_mdu_info <- renderText({
                t <- cem_reactive()
                #t <- ifelse(test = input$cem_approach == "Benevolent" , benevolent_cem_reactive() , aggressive_cem_reactive() ) 
                t<- unfolding(delta = round((1-t),2),ndim = 2)
                t$stress
        })
        
        #download the graph
        output$download_cem_unfolding <- downloadHandler(
                filename = "CEM_MDU.png",
                content = function(file) {
                        #png(file)
                        #print(cem_unfolding_plot())
                        #dev.off()
                        ggsave(file, plot = cem_unfolding_plot(), device = "png", dpi = 450)
                         }
                ) 
        
#End of CEM MDU
        
# Adler Co-plot       
#####
        
        
})