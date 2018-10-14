fluidPage(
  titlePanel("Vaccine Prediction"),
  fluidRow(
    column(3, wellPanel(
      textInput("mhc", "MHC:"),
      textInput("seq","Sequence:"),
      submitButton("Submit")
    )),
    column(6,
           verbatimTextOutput("pep")
    ),
    column(7,verbatimTextOutput("acc")),
    column(8,verbatimTextOutput("accu"))
  )
)
