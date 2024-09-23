#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


# load in packages
library(readxl) # for reading
library(tidyverse)
library(openxlsx) # for writing
#library(parallel)
#library(ggsignif)
#library(ggpattern)
#library(colorspace)
library(ggrepel)
library(shiny)
library(readxl)
#library(vroom)
library(shinyFeedback)
#library(gghighlight)
#library(shinylive)
#library(httpuv)


#setwd("~/Library/CloudStorage/Box-Box/Scripts/Hope_Scripts/Predictions")

#rm(list=ls())

source("volcano.R")
source("combo_trends.R")

# import data
class_colors <- read.csv("class_colors_18Sep2024.csv")
drugs <- read.csv("drug_classes_18Sep2024.csv") %>%
  mutate(fullName = paste0(name, " (", abbreviation, ")"))
input <- "2024-08-23_predictions_Aug2024Data_GRMax_BPaL.xlsx"
volcanoDataset <- read_xlsx("volcano_3_way_august_2024.xlsx") %>%
  mutate(n = as.numeric(n)) %>%
  mutate(n_wo = as.numeric(n_wo)) %>%  
  mutate(prob_with = as.numeric(prob_with)) %>%
  mutate(prob_wo = as.numeric(prob_wo)) %>%
  mutate(prob_diff = as.numeric(prob_diff)) %>%
  mutate(P_wilcox = as.numeric(P_wilcox)) %>%
  mutate(P_wilcox_log10 = as.numeric(P_wilcox_log10)) 

#define lists of interest
drugList <- drugs$abbreviation

classList <- sort(unique(drugs$class))

seriesList <-sort(unique(drugs$series))

#define thresholds
fc_thresh <- 1.5
sig_thresh <- 0.05

#format app dataset
appData <- volcanoDataset %>%
  mutate(signif_pair = cut(P_wilcox, breaks = c(0, sig_thresh, Inf))) %>%
  mutate(signif_pair = case_match(signif_pair,
                                  paste0("(0,", sig_thresh, "]") ~ "Signif",
                                  paste0("(", sig_thresh, ",Inf]") ~ "Not signif")) %>%
  mutate(label = case_when(
    signif_pair == "Signif" & fold_change_log2 >= log2(fc_thresh) ~ "Outperforms",
    signif_pair == "Signif" & fold_change_log2 < log2(1/fc_thresh) ~ "Underperforms",
    .default = "Not significant"
  )) #%>%
# drop_na()


#Add in individual drug, class, and series columns 
for (i in c(1:3)) {
  appData <- appData %>%
    mutate("Drug{i}":= sapply(strsplit(appData$subcombo, split = "\\+"), "[[", i))
  
  appData <- appData %>%
    mutate("Class{i}":= drugs$class[match(sapply(strsplit(appData$subcombo, split = "\\+"), "[[", i), drugs$abbreviation)])
  
  appData <- appData %>%
    mutate("Series{i}":= drugs$series[match(sapply(strsplit(appData$subcombo, split = "\\+"), "[[", i), drugs$abbreviation)])
}


#merge appData with predictions file 
predictions <- read_xlsx(input) %>%
  filter(NumbDrugs == 3)

# App itself!

#Current features
#1) Highlights points containing one drug of interest (DONE 10Jan2024)
#2) Generates a table including the combination, axes' values, drug label for 
#   each point  a user clicks on (DONE 10Jan2024)
#3) Add in full drug list (DONE 12Jan2024 with bugs (see below))
#4) Add in ability to view two drugs at once (DONE 12Jan2024)
#5) Add in drug class and series to the overall data table (DONE 12Jan2024)
#6) Add option to sort drugs visible by class (DONE 12Jan2024)
#7) Fixed string matching issue (i.e., M35 results showing up for M3) (DONE 22Jan2024)
#8) Fixed issue where not all drugs in a category appeared as an option (DONE 22Jan2024)

#Issues to fix
#1) Create more intuitive labels
#2) Add option for how many drugs you want to view at once
#3) Add tab for viewing by class instead of by drug

#Features to add
#1) Add P>BPaL to app Data table and incorporate into table
#2) Output downloadable table with P>BPaL and performance label
#3) Incorporate 3-way backbones when possible


ui <- fluidPage(
  
  #Create title ----------------------------------------------------------------
  
  titlePanel("TB DiaMOND volcano plot"),
  
  fluidRow(
    
    column(3, 
           
           #Create input buttons for drug 1 -------------------------------------------
           
           radioButtons("class1", label = "First drug target of interest:", choices = classList),
           
           selectInput("drug1", label = "First drug of interest:", choices = NULL),
           
           #Create input buttons for drug 2 -------------------------------------------
           
           radioButtons("class2", label = "Second drug target of interest:", choices = classList),
           
           selectInput("drug2", label = "Second drug of interest:", choices = NULL),
           
    ),
    
    column(9, 
           
           #Create tables and plot ----------------------------------------------
           
           plotOutput("volcano", click = "plot_click"),
           
           #abbreviation table
           
           tableOutput("data"),
           
           tableOutput("abbs"),

           
    )
  ),
  
  
  #Set up shinyFeedback ------------------------------------------------------
  
  shinyFeedback::useShinyFeedback(),
  
  
)


server <- function(input, output, session) {
  
  #Observe series selection ----------------------------------------------------
  
  freezeReactiveValue(input, "drug1")
  
  observeEvent(input$class1, {
    choices1 <- drugs %>%
      filter(class == input$class1) %>%
      select(abbreviation)
    updateSelectInput(inputId = "drug1", choices = choices1)
  })
  
  freezeReactiveValue(input, "drug2")
  
  observeEvent(input$class2, {
    choices2 <- drugs %>%
      filter(class == input$class2) %>%
      select(abbreviation)    
    updateSelectInput(inputId = "drug2", choices = choices2)
  })
  
  #set up plot reactives -------------------------------------------------------
  
  drugInput <- reactive({
    appData %>%
      mutate(contains_drugInt1 = input$drug1 == Drug1 | input$drug1 == Drug2 | input$drug1 == Drug3) %>%
      mutate(contains_drugInt2 = input$drug2 == Drug1 | input$drug2 == Drug2 | input$drug2 == Drug3)
  })
  
  clickResult <- reactive({
    appData %>%
      select(subcombo, fold_change_log2, P_wilcox_log10, label)
  })
  
  
  #Output volcano plot with drug of interest -----------------------------------
  
  output$volcano <- renderPlot({
    
    ggplot(drugInput(),
           aes(
             x = fold_change_log2,
             y = P_wilcox_log10,
             color = interaction(contains_drugInt1, contains_drugInt2))) +
      geom_point(alpha = 0.3) +
      theme_classic() +
      geom_vline(xintercept = c(log2(1/fc_thresh), log2(fc_thresh)), col = "brown4", linewidth = 0.2) +
      geom_hline(yintercept = -log10(sig_thresh), col = "brown4", linewidth = 0.2) +
      scale_color_manual(labels = c("Neither drug of interest", input$drug1,
                                    input$drug2, paste0(input$drug1, " and ", input$drug2)),
                         values = c("lightgrey", "blue", "red", "purple"),
                         name = "Combination contains:") +
      xlab("log2(fold change P>BPaL +/- combo in 4 way)") +
      ylab("-log10(P value)") + 
      theme(legend.text = element_text(size=13),
            legend.title = element_text(size=13))
  })
  
  
  #Output table of combinations when clicked -----------------------------------
  
  output$data <- renderTable({
    req(input$plot_click)
    nearPoints(clickResult(), input$plot_click)
  })
  
  freezeReactiveValue(input, "plot_click")
  
  #Output relavent table of abbreviations --------------------------------------
  
  output$abbs <- renderTable({
    drugs %>%
      filter(class == input$class1 | class == input$class2) %>%
      select(name, abbreviation, class)
  })
  
  
  #Stop app when closed (WILL NEED TO BE CHANGED WHEN APP IS DEPLOYED)----------
  
  session$onSessionEnded(function() {
    stopApp()
  })
}


shinyApp(ui, server)
