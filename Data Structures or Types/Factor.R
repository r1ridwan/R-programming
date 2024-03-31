## Factor
# factor is a data type used to represent categorical data
smooking_statrus<- c("yes", "no", "no", "yes")
factor(smooking_statrus)
factor(c("yes", "no", "no", "yes"))


# Create a factor with unordered levels
gender <- factor(c("Male", "Female", "Male", "Male", "Female"))
# Display the factor
gender
# print the levels
levels(gender)


# factor with custom levels
gender2 <- factor(c("Male", "Female", "Male", "Male", "Female"),
                  ordered = T,
                  labels = c("M", "F"))


# Ordered factors
age_group <- factor(c("Young", "Middle-aged", "Senior", "Young", "Senior"), 
                    ordered = TRUE,
                    levels = c("Young", "Middle-aged", "Senior"))






