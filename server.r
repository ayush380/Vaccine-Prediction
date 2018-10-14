function(input, output) {
  
  # output$plot1 <- renderPlot({
  #   hist(rnorm(input$n))
  # })
  # output$wait<- renderText({paste("Please Wait")})

  output$pep <- renderText({
    source('./src/driver.R')
    #setwd("C:/Users/Ayush Srivastava/Desktop/Minor")
    req(input$mhc,input$seq,cancelOutput = TRUE)
    
    #paste(" ")
    paste(input$mhc,input$seq)
    lala<-driver(input$mhc,input$seq)
    c<-paste(" ",lala)
    
    #c<-paste(c,"/n",'Random Forest,Smm, Prbability',"/n",lala$accuracies)
  })
  
  # output$seq <- renderText({
  #   paste("Seq text is:", input$seq)
  # })
  
  output$accu<-renderText({
    req(input$mhc,input$seq)
    lala<-Acc(input$mhc,input$seq)
    paste(lala)
    
  }
  )
  output$acc<-renderText({
    paste("Accuracies\n","Random Forest","SMM","Probability")
  })
  
  
  
}