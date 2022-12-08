
rm(list=ls())
set.seed( 2011 )
MWade = 135.13 # molecular weights for adenine in g/mol
MWcre = 113.12 #molecular weights for creatinine in g/mol

## create the dataset
# GRPA1 group
GRP1<-tibble::tribble(
  ~Variable,   ~Level,        ~Number,         ~`HR.(univariable)`,       ~`HR.(multivariable)`,
  "Sex", "Female", "2204 (100.0)",                          NA,                          NA,
  NA,   "Male", "2318 (100.0)", "1.13 (0.91-1.40, p=0.265)", "1.13 (0.91-1.40, p=0.276)",
  "Score",      "1", "2401 (100.0)",                          NA,                          NA,
  NA,    "1-2", "1637 (100.0)", "1.49 (1.19-1.87, p=0.001)", "1.15 (0.90-1.47, p=0.250)",
  NA,    "3-4",  "412 (100.0)", "1.71 (1.14-2.56, p=0.010)", "1.09 (0.71-1.67, p=0.710)",
  NA,    ">=5",   "42 (100.0)", "1.67 (0.53-5.21, p=0.381)", "0.96 (0.30-3.05, p=0.943)",
  "Treatment",      "A", "1572 (100.0)",                          NA,                          NA,
  NA,      "B", "2951 (100.0)", "1.74 (1.26-2.40, p=0.001)", "1.53 (1.09-2.13, p=0.013)"
)

# GRP2 group
GRP2<-tibble::tribble(
  ~Variable,   ~Level,        ~Number,          ~`HR.(univariable)`,          ~`HR.(multivariable)`,
  "Sex", "Female", "2204 (100.0)",                           NA,                           NA,
  NA,   "Male", "2318 (100.0)",  "1.70 (1.36-2.13, p<0.001)",  "1.62 (1.28-2.04, p<0.001)",
  "Score",      "1", "2401 (100.0)",                           NA,                           NA,
  NA,    "1-2", "1637 (100.0)",  "2.76 (1.21-6.29, p=0.016)",  "2.69 (1.18-6.13, p=0.019)",
  NA,    "3-4",  "412 (100.0)", "5.11 (2.26-11.58, p<0.001)", "4.46 (1.95-10.23, p<0.001)",
  NA,    ">=5",   "42 (100.0)", "5.05 (2.19-11.64, p<0.001)",  "4.08 (1.73-9.59, p=0.001)",
  "Treatment",      "A", "1572 (100.0)",                           NA,                           NA,
  NA,      "B", "2951 (100.0)",  "1.48 (1.16-1.88, p=0.001)",  "1.23 (0.95-1.59, p=0.114)"
)




### creating function for desired column for the forestplot. 

library(dplyr)

wrangler <- function(data){
  grp <- as.character(match.call()$data)
  data %>%
    tidyr::fill(Variable) %>% 
    mutate(Variable = paste(Variable, Level), 
           Number = as.numeric(gsub("^(\\d+).*$", "\\1", Number)),
           univariable_HR =  as.numeric(gsub("^((\\d+|\\.)+).*$", "\\1", `HR.(univariable)`)),
           univariable_lower = as.numeric(gsub("^.+? \\((.+?)-.*$", "\\1", `HR.(univariable)`)),
           univariable_upper = as.numeric(gsub("^.+?-(.+?),.*$", "\\1", `HR.(univariable)`)),
           univariable_p = gsub("^.+?p=*(.+?)\\).*$", "\\1", `HR.(univariable)`),
           multivariable_HR =  as.numeric(gsub("^((\\d+|\\.)+).*$", "\\1", `HR.(multivariable)`)),
           multivariable_lower = as.numeric(gsub("^.+? \\((.+?)-.*$", "\\1", `HR.(multivariable)`)),
           multivariable_upper = as.numeric(gsub("^.+?-(.+?),.*$", "\\1", `HR.(multivariable)`)),
           multivariable_p = gsub("^.+?p=*(.+?)\\).*$", "\\1", `HR.(multivariable)`),
           group = grp) %>% 
    filter(!is.na(univariable_HR)) %>%
    select(-Level, -`HR.(multivariable)`, - `HR.(univariable)`) %>%
    tidyr::pivot_longer(cols = -(c(1:2, 11)), names_sep = "_", names_to = c("type", ".value"))
}

### Passing the two growp into the function

wrangler(GRP1)
wrangler(GRP2)

# Rowbind two dataset
df <- rbind(wrangler(GRP1), wrangler(GRP2))

###  Head of the main dataset

head(df)

getwd()

### Saving the main dataset of our directory

write.csv(df,"C:/Users/USER/Desktop/Rabiul Sir/New work/data.csv")

### Visualization ggplot2

library(ggplot2)

### Simple point forest plot with ggplot2####

 ggplot(df, aes(HR, Variable)) +
  geom_pointrange(aes(xmin = lower, xmax = upper, colour = type),
                  position = position_dodge(width = 0.5)) +
  facet_grid(group~., switch = "y")  + 
 #geom_errorbar(aes(xmin=lower, xmax=upper,col=type),width=0.5,cex=1)+ 

  geom_vline(xintercept = 0, linetype = 2) +
 #theme_bw() +
  theme(
    strip.placement = "outside",
        strip.text= element_text(angle = 180),
        strip.background = element_blank(),
        panel.spacing = unit(0, "mm"))

### Error bar forest plot with ggplot2####

ggplot(df, aes(HR, Variable)) +
  geom_pointrange(aes(xmin = lower, xmax = upper, colour = type),
                  position = position_dodge(width = 0.5)) +
  facet_grid(group~., switch = "y")  + 
  geom_errorbar(aes(xmin=lower, xmax=upper,col=type),width=0.5,cex=0.5)+ 
  
  geom_vline(xintercept = 0, linetype = 2) +
  #theme_bw() +
  theme(
    strip.placement = "outside",
    strip.text= element_text(angle = 180),
    strip.background = element_blank(),
    panel.spacing = unit(0, "mm"))









  











