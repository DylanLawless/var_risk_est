# https://www.iaincollings.com/probability-and-random-variables

library(ggplot2)
library(gridExtra)

mu <- 0
Changing <- 1  # 1 for changing sig, 2 for changing theta
if (Changing == 1) {
  theta <- 1
  ScenarioN <- 5
} else {
  thetaStart <- 1
  sig <- 2
  ScenarioN <- 2
  ThetaDiffFactor <- 1
}
LowVal <- -5
HighVal <- 10
Res <- 100
n <- seq(LowVal, HighVal, by = 1/Res)

data <- data.frame()
for (Scenario in 1:ScenarioN) {
  if (Changing == 1) {
    sig <- (Scenario + 3) / 20
    currentTheta <- theta
  } else {
    currentTheta <- thetaStart + (Scenario - 1) * ThetaDiffFactor
  }
  f <- 1 / (sqrt(2 * pi) * sig) * exp( - (n - mu - currentTheta)^2 / (2 * sig^2))
  Lnf <- log(f)
  DerLnf <- (n - mu - currentTheta) / (sig^2)
  SqDerLnf <- DerLnf^2
  fSqDerLnf <- f * SqDerLnf
  tmp <- data.frame(
    Scenario = factor(Scenario),
    n = n,
    f = f,
    Lnf = Lnf,
    DerLnf = DerLnf,
    SqDerLnf = SqDerLnf,
    fSqDerLnf = fSqDerLnf
  )
  data <- rbind(data, tmp)
}

p1 <- ggplot(data, aes(x = n, y = f, colour = Scenario)) +
  geom_line() +
  labs(title = "PDF", x = "n", y = "f") +
  theme_minimal()

p2 <- ggplot(data, aes(x = n, y = Lnf, colour = Scenario)) +
  geom_line() +
  labs(title = "log(PDF)", x = "n", y = "log(f)") +
  theme_minimal()

p3 <- ggplot(data, aes(x = n, y = DerLnf, colour = Scenario)) +
  geom_line() +
  labs(title = "d/dtheta", x = "n", y = "d/dtheta") +
  theme_minimal()

p4 <- ggplot(data, aes(x = n, y = SqDerLnf, colour = Scenario)) +
  geom_line() +
  labs(title = "(d/dtheta)^2", x = "n", y = "(d/dtheta)^2") +
  theme_minimal()

p5 <- ggplot(data, aes(x = n, y = fSqDerLnf, colour = Scenario)) +
  geom_line() +
  labs(title = "(d/dtheta)^2 PDF", x = "n", y = "(d/dtheta)^2 * PDF") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, ncol = 1)









# data ----


source("2025_iuis_iei_table.R")




set.seed(123)  # for reproducibility
n <- 100  # number of rows

df <- data.frame(
  `Major category` = sample(c("CID", "PID", "Other"), n, replace = TRUE),
  Subcategory = sample(paste("Sub", 1:5), n, replace = TRUE),
  Disease = sample(c("Disease A", "Disease B", "Disease C"), n, replace = TRUE),
  `Genetic defect` = sample(paste("Gene", LETTERS), n, replace = TRUE),
  Inheritance = sample(c("AR", "AD", "X-linked"), n, replace = TRUE),
  `OMIM ID` = sample(100000:200000, n, replace = TRUE),
  Uniprot = sample(paste0("P", sprintf("%05d", sample(1:99999, n, replace = TRUE))), n, replace = TRUE),
  VariantCounts = paste0(sample(10:50, n, replace = TRUE), " / ", sample(300:600, n, replace = TRUE)),
  `T cell` = sample(1:5, n, replace = TRUE),
  `B cell` = sample(1:5, n, replace = TRUE),
  Ig = sample(1:5, n, replace = TRUE),
  Neutrophil = sample(1:5, n, replace = TRUE),
  `T cell details` = sample(c("detail A", "detail B", "detail C", "detail D"), n, replace = TRUE),
  `B cell details` = sample(c("detail A", "detail B", "detail C", "detail D"), n, replace = TRUE),
  `Ig details` = sample(c("detail A", "detail B", "detail C", "detail D"), n, replace = TRUE),
  `Neutrophil details` = sample(c("detail A", "detail B", "detail C", "detail D"), n, replace = TRUE),
  `Associated features` = sample(c("feature 1", "feature 2", "feature 3"), n, replace = TRUE),
  ICD9 = sample(100:999, n, replace = TRUE),
  ICD10 = sample(paste0("A", sample(100:999, n, replace = TRUE)), n, replace = TRUE),
  `HPO IDs` = sample(1000:9999, n, replace = TRUE),
  `HPO term` = sample(c("term1", "term2", "term3", "term4"), n, replace = TRUE),
  `Other affected cells` = sample(c("cell1", "cell2", "cell3"), n, replace = TRUE),
  `Inheritance detail` = sample(c("detail1", "detail2", "detail3"), n, replace = TRUE),
  `GOF/DN details` = sample(c("detail1", "detail2", "detail3"), n, replace = TRUE),
  Major_category_original = sample(c("CID", "PID", "Other"), n, replace = TRUE),
  AlphaFold_URL = paste0("http://alphafold.ebi.ac.uk/", sample(letters, n, replace = TRUE)),
  Other = sample(c("misc1", "misc2", "misc3"), n, replace = TRUE),
  Pathogenic = sample(c("Yes", "No"), n, replace = TRUE),
  stringsAsFactors = FALSE
)

head(df)






# Set up fake data and simulate variant frequency–based outcomes for two features:
set.seed(123)
n <- 200

# Use the previously defined fake dataset (df) or create a new one for this simulation:
df <- data.frame(
  `T cell` = sample(1:5, n, replace = TRUE),
  `T cell details` = sample(c("detail A", "detail B", "detail C", "detail D"), n, replace = TRUE)
)

# Convert column names to valid R names
names(df) <- make.names(names(df))

# Create numeric representations:
# Structured feature ("T cell") is already numeric.
# For unstructured "T cell details", assign numeric codes.
df$T_cell_details_numeric <- as.numeric(factor(df$T.cell.details, 
                                               levels = c("detail A", "detail B", "detail C", "detail D")))

# Simulate binary outcomes (Pathogenic) via logistic models.
# For the structured feature "T cell" we assume a stronger effect (beta = 1)
# so that the outcome depends more clearly on this feature.
logit_prob1 <- 0 + 1 * df$T.cell  # intercept 0, beta 1
p1 <- exp(logit_prob1) / (1 + exp(logit_prob1))
df$Pathogenic_bin <- rbinom(n, 1, p1)

# For the unstructured feature ("T cell details") we assume a weaker effect (beta = 0.5)
logit_prob2 <- 0 + 0.5 * df$T_cell_details_numeric
p2 <- exp(logit_prob2) / (1 + exp(logit_prob2))
df$Pathogenic_bin2 <- rbinom(n, 1, p2)

# Fit logistic regression models for each feature:
model_struct <- glm(Pathogenic_bin ~ T.cell, data = df, family = binomial(link = "logit"))
model_unstruct <- glm(Pathogenic_bin2 ~ T_cell_details_numeric, data = df, family = binomial(link = "logit"))

# Display summaries:
summary(model_struct)
summary(model_unstruct)

# Compute the observed Fisher information for the beta (slope) parameter.
# In logistic regression, the variance of the MLE for beta is approximately the inverse of the observed information.
# Thus, observed Fisher information can be estimated as:
# I(beta) = 1 / Var(beta.hat)
fisher_info_struct <- 1 / (summary(model_struct)$coefficients[2, 2]^2)
fisher_info_unstruct <- 1 / (summary(model_unstruct)$coefficients[2, 2]^2)

cat("Fisher information for 'T cell' feature (structured):", fisher_info_struct, "\n")
cat("Fisher information for 'T cell details' feature (unstructured):", fisher_info_unstruct, "\n")









library(ggplot2)
library(gridExtra)

# Define a function to compute the likelihood and related quantities for logistic regression
# with a fixed intercept (here, set to 0 for simplicity)
compute_likelihood <- function(beta, x, y) {
  # For each beta in the grid, compute the log-likelihood:
  # logL(beta) = sum_i [ y_i*(beta * x_i) - log(1 + exp(beta * x_i)) ]
  logL <- sapply(beta, function(b) sum( y * (b * x) - log(1 + exp(b * x)) ) )
  L <- exp(logL)
  
  # Derivative of the log-likelihood with respect to beta:
  # d/d_beta logL(beta) = sum_i [ x_i * (y_i - p_i) ]
  # where p_i = exp(b * x_i) / (1 + exp(b * x_i))
  dlogL <- sapply(beta, function(b) {
    p <- exp(b * x) / (1 + exp(b * x))
    sum( x * (y - p) )
  })
  
  SqDerLnf <- dlogL^2
  fSqDerLnf <- L * SqDerLnf
  
  data.frame(beta = beta, L = L, logL = logL,
             dlogL = dlogL, SqDerLnf = SqDerLnf, fSqDerLnf = fSqDerLnf)
}

# Create a grid for the parameter beta
beta_grid <- seq(-3, 3, length.out = 200)

# For demonstration, we use the earlier simulated df.
# Our features are 'T.cell' (structured) and 'T_cell_details_numeric' (unstructured)
# and the outcomes are simulated as Pathogenic_bin and Pathogenic_bin2 respectively.
# (Assuming df has been defined as in the previous code block.)

# Extract predictor (x) and outcome (y) for each feature:
x_struct <- df$T.cell
y_struct <- df$Pathogenic_bin

x_unstruct <- df$T_cell_details_numeric
y_unstruct <- df$Pathogenic_bin2

# Compute the likelihood functions for both features:
likelihood_struct <- compute_likelihood(beta_grid, x_struct, y_struct)
likelihood_struct$Feature <- "T cell"

likelihood_unstruct <- compute_likelihood(beta_grid, x_unstruct, y_unstruct)
likelihood_unstruct$Feature <- "T cell details"

# Combine the data for plotting
likelihood_data <- rbind(likelihood_struct, likelihood_unstruct)

# Create analogous plots to the previous example:

p1 <- ggplot(likelihood_data, aes(x = beta, y = L, colour = Feature)) +
  geom_line() +
  labs(title = "Likelihood", x = expression(beta), y = expression(L(beta))) +
  theme_minimal()

p2 <- ggplot(likelihood_data, aes(x = beta, y = logL, colour = Feature)) +
  geom_line() +
  labs(title = "Log-Likelihood", x = expression(beta), y = expression(log~L(beta))) +
  theme_minimal()

p3 <- ggplot(likelihood_data, aes(x = beta, y = dlogL, colour = Feature)) +
  geom_line() +
  labs(title = "Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression(frac(d~log~L, d~beta))) +
  theme_minimal()

p4 <- ggplot(likelihood_data, aes(x = beta, y = SqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = "Squared Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression((frac(d~log~L, d~beta))^2)) +
  theme_minimal()

p5 <- ggplot(likelihood_data, aes(x = beta, y = fSqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2),
       x = expression(beta), 
       y = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2)) +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, ncol = 1)







Assuming \( p_{\text{pathogenic}} \) is the probability for a given individual to have a pathogenic variant in that gene, then the probability of observing at least one patient in a population of size \( N \) is

\[
  P(\text{at least one}) = 1 - (1 - p_{\text{pathogenic}})^N.
  \]

For example, if \( p_{\text{pathogenic}} = 1 - (1 - 0.005)^{10} \) for one gene and \( N = 1,000,000 \), then

\[
  P(\text{at least one}) = 1 - (1 - p_{\text{pathogenic}})^{1,000,000}.
  \]

Below is an R code snippet that calculates this:
  

# Parameters for a single gene:
n_possible <- 10000       # total possible SNVs per gene
n_pathogenic <- 1     # known pathogenic SNVs
frequency <- 0.005      # frequency of each pathogenic SNV (0.5%)

# Probability for one individual:
p_pathogenic <- 1 - (1 - frequency)^n_pathogenic

# Population size:
population_size <- 500000

# Probability of observing at least one patient in the population with the variant:
prob_at_least_one <- 1 - (1 - p_pathogenic)^population_size

cat("Probability of observing at least one patient with the variant:", 
    round(prob_at_least_one, 4), "\n")


This calculation gives you the chance that in a given population, at least one patient has a pathogenic variant in that gene.




















For an autosomal dominant (AD) gene where one copy of the pathogenic variant is sufficient, the calculation

```r
p_pathogenic <- 1 - (1 - frequency)^n_pathogenic
```

with `n_pathogenic <- 1` and `frequency <- 0.005` correctly gives

\[
  p_{AD} = 1 - (1-0.005)^1 = 0.005.
  \]

However, for an autosomal recessive (AR) gene you need two copies of a pathogenic variant for the disease to manifest. Under Hardy–Weinberg assumptions (and assuming the variant is rare) the probability that a given individual has two copies of the same pathogenic allele is approximately

\[
  p_{AR} \approx (\text{frequency})^2.
  \]

If you have only one variant that can cause disease, then with `frequency <- 0.005` the probability becomes

```r
p_AR <- frequency^2
```

which is

\[
  0.005^2 = 0.000025.
  \]

If there are multiple independent pathogenic variants (say \( n \) variants) in the gene that can each cause disease via compound heterozygosity or homozygosity, one common approximation is to first estimate the probability that a randomly chosen allele is pathogenic as

\[
  p_{\text{allele}} = n_{\text{pathogenic}} \times \text{frequency}
  \]

(assuming the variants are rare and non-overlapping). Then, under Hardy–Weinberg, the probability that an individual carries two pathogenic alleles is approximately

\[
  p_{AR} \approx (p_{\text{allele}})^2 = (n_{\text{pathogenic}} \times \text{frequency})^2.
  \]

For example, if \( n_{\text{pathogenic}} = 10 \) and `frequency <- 0.005`, then

```r
n_pathogenic <- 10
frequency <- 0.005
p_AR <- (n_pathogenic * frequency)^2
p_AR
```

This gives

\[
  (10 \times 0.005)^2 = (0.05)^2 = 0.0025.
  \]

Thus, in your analysis you need to choose the appropriate calculation based on the inheritance pattern. For AD, one pathogenic variant is enough, so the probability is \(0.005\) per allele; for AR, the chance is the square of the allele pathogenicity probability, which might be \( (0.005)^2 \) for one variant or \((10 \times 0.005)^2\) if considering 10 independent variants.









# new ----

library(ggplot2)
library(gridExtra)

# --- Data Simulation with Inheritance Patterns ---

set.seed(123)
n <- 200

# Function to generate a random string of letters of random length between 3 and 10:
random_string <- function() {
  len <- sample(3:10, 1)
  paste(sample(letters, len, replace = TRUE), collapse = "")
}

# Simulate features and inheritance patterns:
df <- data.frame(
  T.cell = sample(1:5, n, replace = TRUE),
  # Generate random text descriptions for T.cell.details:
  T.cell.details = replicate(n, random_string()),
  Inheritance = sample(c("AD", "AR", "X-linked"), n, replace = TRUE),
  stringsAsFactors = FALSE
)

# Convert column names to valid R names for later use
names(df) <- make.names(names(df))

# Create a numeric representation for the unstructured feature by converting the text to a factor
# Note: This numeric representation is solely for simulation purposes and does not capture the semantics of text.
df$T_cell_details_numeric <- as.numeric(factor(df$T.cell.details))




# text ----

library(text2vec)
library(irlba)

# Create an iterator over the T.cell.details column
it <- itoken(df$T.cell.details, 
             tokenizer = word_tokenizer, 
             progressbar = FALSE)

# Build a vocabulary and vectorizer
vocab <- create_vocabulary(it)
vectorizer <- vocab_vectorizer(vocab)

# Create a Document-Term Matrix (DTM)
dtm <- create_dtm(it, vectorizer)

# Apply TF-IDF weighting
tfidf_transformer <- TfIdf$new()
dtm_tfidf <- tfidf_transformer$fit_transform(dtm)

# Perform SVD for dimensionality reduction (e.g. reduce to 1 dimension)
# This produces a more semantically-informed numeric representation.
svd_result <- irlba(dtm_tfidf, nv = 1)
df$T_cell_details_embedding <- svd_result$u[, 1]

# Check the new embedding:
head(df$T_cell_details_embedding)




# --- Calculate Baseline Pathogenic Probability Based on Inheritance ---

frequency <- 0.005  # frequency for a single pathogenic SNV

df$p_base <- ifelse(df$Inheritance %in% c("AD", "X.linked"),
                    frequency,      # AD / X-linked: one copy is enough
                    frequency^2)    # AR: need two copies

df$logit_base <- log(df$p_base / (1 - df$p_base))

# --- Simulate Binary Outcomes Using Logistic Models with Feature Effects ---

# For the structured feature "T.cell", assume a stronger effect (beta = 1)
logit_prob1 <- df$logit_base + 1 * df$T.cell
p1 <- exp(logit_prob1) / (1 + exp(logit_prob1))
df$Pathogenic_bin <- rbinom(n, 1, p1)

# For the unstructured feature "T.cell.details" (using its numeric representation), assume a weaker effect (beta = 0.5)
logit_prob2 <- df$logit_base + 0.5 * df$T_cell_details_numeric
p2 <- exp(logit_prob2) / (1 + exp(logit_prob2))
df$Pathogenic_bin2 <- rbinom(n, 1, p2)

# --- Fit Logistic Regression Models ---

model_struct <- glm(Pathogenic_bin ~ T.cell, data = df, family = binomial(link = "logit"))
model_unstruct <- glm(Pathogenic_bin2 ~ T_cell_details_numeric, data = df, family = binomial(link = "logit"))

summary(model_struct)
summary(model_unstruct)

# Compute observed Fisher information for the slope parameter (beta)
fisher_info_struct <- 1 / (summary(model_struct)$coefficients[2, 2]^2)
fisher_info_unstruct <- 1 / (summary(model_unstruct)$coefficients[2, 2]^2)

cat("Fisher information for 'T.cell' (structured):", fisher_info_struct, "\n")
cat("Fisher information for 'T.cell.details' (unstructured):", fisher_info_unstruct, "\n")

# --- Compute Likelihood Functions and Related Quantities ---

compute_likelihood <- function(beta, x, y) {
  logL <- sapply(beta, function(b) sum( y * (b * x) - log(1 + exp(b * x)) ) )
  L <- exp(logL)
  dlogL <- sapply(beta, function(b) {
    p <- exp(b * x) / (1 + exp(b * x))
    sum( x * (y - p) )
  })
  SqDerLnf <- dlogL^2
  fSqDerLnf <- L * SqDerLnf
  data.frame(beta = beta, L = L, logL = logL,
             dlogL = dlogL, SqDerLnf = SqDerLnf, fSqDerLnf = fSqDerLnf)
}

beta_grid <- seq(-3, 3, length.out = 200)

likelihood_struct <- compute_likelihood(beta_grid, df$T.cell, df$Pathogenic_bin)
likelihood_struct$Feature <- "T.cell"

likelihood_unstruct <- compute_likelihood(beta_grid, df$T_cell_details_numeric, df$Pathogenic_bin2)
likelihood_unstruct$Feature <- "T.cell.details"

likelihood_data <- rbind(likelihood_struct, likelihood_unstruct)

# --- Create Likelihood Plots ---

p1 <- ggplot(likelihood_data, aes(x = beta, y = L, colour = Feature)) +
  geom_line() +
  labs(title = "Likelihood", x = expression(beta), y = expression(L(beta))) +
  theme_minimal()

p2 <- ggplot(likelihood_data, aes(x = beta, y = logL, colour = Feature)) +
  geom_line() +
  labs(title = "Log-Likelihood", x = expression(beta), y = expression(log~L(beta))) +
  theme_minimal()

p3 <- ggplot(likelihood_data, aes(x = beta, y = dlogL, colour = Feature)) +
  geom_line() +
  labs(title = "Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression(frac(d~log~L, d~beta))) +
  theme_minimal()

p4 <- ggplot(likelihood_data, aes(x = beta, y = SqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = "Squared Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression((frac(d~log~L, d~beta))^2)) +
  theme_minimal()

p5 <- ggplot(likelihood_data, aes(x = beta, y = fSqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2),
       x = expression(beta), 
       y = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2)) +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, ncol = 1)

# --- Population-Level Probability Example ---

n_possible <- 10000  # total possible SNVs per gene
inheritance_type <- "AD"  # choose an inheritance type for demonstration

if (inheritance_type %in% c("AD", "X-linked")) {
  p_pathogenic <- frequency
} else if (inheritance_type == "AR") {
  p_pathogenic <- frequency^2
}

population_size <- 500000
prob_at_least_one <- 1 - (1 - p_pathogenic)^population_size

cat("For inheritance type", inheritance_type, "\n")
cat("Probability of an individual carrying the pathogenic variant:", p_pathogenic, "\n")
cat("Probability of observing at least one patient in the population:", 
    round(prob_at_least_one, 4), "\n")






















# last ----





library(ggplot2)
library(gridExtra)
library(text2vec)
library(irlba)

set.seed(123)
n <- 200

# Function to generate a random string of letters of random length between 3 and 10:
random_string <- function() {
  len <- sample(3:10, 1)
  paste(sample(letters, len, replace = TRUE), collapse = "")
}

# --- Data Simulation with Inheritance Patterns ---

df <- data.frame(
  T.cell = sample(1:5, n, replace = TRUE),
  # Generate random text descriptions for T.cell.details:
  T.cell.details = replicate(n, random_string()),
  Inheritance = sample(c("AD", "AR", "X-linked"), n, replace = TRUE),
  stringsAsFactors = FALSE
)

# Convert column names to valid R names for later use
names(df) <- make.names(names(df))

# Create a numeric representation for the unstructured feature using factor conversion:
df$T_cell_details_numeric <- as.numeric(factor(df$T.cell.details))

# --- Generate a Text Embedding for T.cell.details using TF-IDF + SVD ---
it <- itoken(df$T.cell.details, 
             tokenizer = word_tokenizer, 
             progressbar = FALSE)
vocab <- create_vocabulary(it)
vectorizer <- vocab_vectorizer(vocab)
dtm <- create_dtm(it, vectorizer)
tfidf_transformer <- TfIdf$new()
dtm_tfidf <- tfidf_transformer$fit_transform(dtm)
svd_result <- irlba(dtm_tfidf, nv = 1)
df$T_cell_details_embedding <- svd_result$u[, 1]

# --- Standardise the Predictors ---
# Standardise the structured feature (T.cell) and the unstructured numeric representations.
df$T_cell_std <- scale(df$T.cell)
df$T_cell_details_numeric_std <- scale(df$T_cell_details_numeric)
df$T_cell_details_embedding_std <- scale(df$T_cell_details_embedding)

# --- Calculate Baseline Pathogenic Probability Based on Inheritance ---

frequency <- 0.005  # frequency for a single pathogenic SNV

df$p_base <- ifelse(df$Inheritance %in% c("AD", "X.linked"),
                    frequency,      # AD / X-linked: one copy is enough
                    frequency^2)    # AR: need two copies
df$logit_base <- log(df$p_base / (1 - df$p_base))

# --- Simulate Binary Outcomes Using Logistic Models with Feature Effects ---
# For the structured feature, assume a stronger effect (beta = 1)
logit_prob1 <- df$logit_base + 1 * df$T_cell_std
p1 <- exp(logit_prob1) / (1 + exp(logit_prob1))
df$Pathogenic_bin <- rbinom(n, 1, p1)

# For the unstructured feature (using the numeric representation), assume a weaker effect (beta = 0.5)
logit_prob2 <- df$logit_base + 0.5 * df$T_cell_details_numeric_std
p2 <- exp(logit_prob2) / (1 + exp(logit_prob2))
df$Pathogenic_bin2 <- rbinom(n, 1, p2)

# Alternatively, if you prefer to use the embedding instead:
# logit_prob2 <- df$logit_base + 0.5 * df$T_cell_details_embedding_std
# p2 <- exp(logit_prob2) / (1 + exp(logit_prob2))
# df$Pathogenic_bin2 <- rbinom(n, 1, p2)

# --- Fit Logistic Regression Models with Standardised Predictors ---

model_struct <- glm(Pathogenic_bin ~ T_cell_std, data = df, family = binomial(link = "logit"))
model_unstruct <- glm(Pathogenic_bin2 ~ T_cell_details_numeric_std, data = df, family = binomial(link = "logit"))

summary(model_struct)
summary(model_unstruct)

# Compute observed Fisher information for the slope parameter (beta)
fisher_info_struct <- 1 / (summary(model_struct)$coefficients[2, 2]^2)
fisher_info_unstruct <- 1 / (summary(model_unstruct)$coefficients[2, 2]^2)

cat("Fisher information for 'T.cell' (structured, standardised):", fisher_info_struct, "\n")
cat("Fisher information for 'T.cell.details' (unstructured, standardised):", fisher_info_unstruct, "\n")

# --- Compute Likelihood Functions and Related Quantities ---

compute_likelihood <- function(beta, x, y) {
  logL <- sapply(beta, function(b) sum( y * (b * x) - log(1 + exp(b * x)) ) )
  L <- exp(logL)
  dlogL <- sapply(beta, function(b) {
    p <- exp(b * x) / (1 + exp(b * x))
    sum( x * (y - p) )
  })
  SqDerLnf <- dlogL^2
  fSqDerLnf <- L * SqDerLnf
  data.frame(beta = beta, L = L, logL = logL,
             dlogL = dlogL, SqDerLnf = SqDerLnf, fSqDerLnf = fSqDerLnf)
}

beta_grid <- seq(-3, 3, length.out = 200)

likelihood_struct <- compute_likelihood(beta_grid, df$T_cell_std, df$Pathogenic_bin)
likelihood_struct$Feature <- "T.cell (std.)"

likelihood_unstruct <- compute_likelihood(beta_grid, df$T_cell_details_numeric_std, df$Pathogenic_bin2)
likelihood_unstruct$Feature <- "T.cell.details (std.)"

likelihood_data <- rbind(likelihood_struct, likelihood_unstruct)

# --- Create Likelihood Plots ---

p1 <- ggplot(likelihood_data, aes(x = beta, y = L, colour = Feature)) +
  geom_line() +
  labs(title = "Likelihood", x = expression(beta), y = expression(L(beta))) +
  theme_minimal()

p2 <- ggplot(likelihood_data, aes(x = beta, y = logL, colour = Feature)) +
  geom_line() +
  labs(title = "Log-Likelihood", x = expression(beta), y = expression(log~L(beta))) +
  theme_minimal()

p3 <- ggplot(likelihood_data, aes(x = beta, y = dlogL, colour = Feature)) +
  geom_line() +
  labs(title = "Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression(frac(d~log~L, d~beta))) +
  theme_minimal()

p4 <- ggplot(likelihood_data, aes(x = beta, y = SqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = "Squared Derivative of Log-Likelihood", 
       x = expression(beta), 
       y = expression((frac(d~log~L, d~beta))^2)) +
  theme_minimal()

p5 <- ggplot(likelihood_data, aes(x = beta, y = fSqDerLnf, colour = Feature)) +
  geom_line() +
  labs(title = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2),
       x = expression(beta), 
       y = expression(L(beta)~"*"~(frac(d~log~L, d~beta))^2)) +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, ncol = 1)

# --- Population-Level Probability Example ---

n_possible <- 10000  # total possible SNVs per gene
inheritance_type <- "AD"  # choose an inheritance type for demonstration

if (inheritance_type %in% c("AD", "X-linked")) {
  p_pathogenic <- frequency
} else if (inheritance_type == "AR") {
  p_pathogenic <- frequency^2
}

population_size <- 500000
prob_at_least_one <- 1 - (1 - p_pathogenic)^population_size

cat("For inheritance type", inheritance_type, "\n")
cat("Probability of an individual carrying the pathogenic variant:", p_pathogenic, "\n")
cat("Probability of observing at least one patient in the population:", 
    round(prob_at_least_one, 4), "\n")







# dbnsfp ----

library(ggplot2)
library(dplyr)

header_line <- readLines("~/Desktop/dbnsfp/data/tnfaip3_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

df <- read.table("~/Desktop/dbnsfp/data/tnfaip3", 
                 sep = "\t",
                 header = FALSE, 
                 stringsAsFactors = FALSE, 
                 fill = TRUE
                 )

colnames(df) <- header_fields

head(df, 1)

hold <- df
df <- hold
# 
# filter ----
df$clinvar_clnsig |> unique()

df <- df |> dplyr::filter(!clinvar_clnsig  == ".")
df|> count(clinvar_clnsig)

ggplot(df, aes(x = clinvar_clnsig)) +
  geom_bar(fill = "skyblue", color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_minimal() +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = "Count of ClinVar Clinical Significance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- NA
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df$gnomAD_genomes_AF

df$`pos(1-based)`

df <- df |> select(`pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig)
df$Inheritance <- "AD"

head(df)

# calc ----
df$clinvar_clnsig |> unique()

# First, filter for known pathogenic variants:
df_path <- subset(df, clinvar_clnsig == "Uncertain_significance")

# Ensure the allele frequency column is numeric:
df_path$gnomAD_genomes_AF <- as.numeric(df_path$gnomAD_genomes_AF)

# Define a fixed population size (you can adjust this as needed)
population_size <- 83702 # births in 2023

# For autosomal dominant (AD) and X-linked, one allele is sufficient,
# so we use the allele frequency directly.
# For autosomal recessive (AR), the probability for an individual is roughly (AF)^2.
df_path$disease_prob <- ifelse(df_path$Inheritance %in% c("AD", "X-linked"),
                               df_path$gnomAD_genomes_AF,
                               df_path$gnomAD_genomes_AF^2)

# Expected number of individuals in the population with the variant:
df_path$expected_cases <- population_size * df_path$disease_prob

# Probability of observing at least one affected individual in the population:
df_path$prob_at_least_one <- 1 - (1 - df_path$disease_prob)^population_size

# View the key results:
df_path[, c("gnomAD_genomes_AF", "Inheritance", "disease_prob", 
            "expected_cases", "prob_at_least_one")]


library(ggplot2)
library(patchwork)

# Filter out rows with NA in allele frequency
df_path_nonNA <- subset(df_path, !is.na(gnomAD_genomes_AF))

# Scatter plot: Expected cases vs Allele frequency
p1 <- ggplot(df_path_nonNA, aes(x = gnomAD_genomes_AF, y = expected_cases, color = Inheritance)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  scale_x_log10() +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Expected Cases",
       title = "Expected Cases vs Allele Frequency") +
  theme_minimal()

# Scatter plot: Probability of at least one case vs Allele frequency
p2 <- ggplot(df_path_nonNA, aes(x = gnomAD_genomes_AF, y = prob_at_least_one, color = Inheritance)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  scale_x_log10() +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Probability of ≥1 Case",
       title = "Probability of at Least One Case vs Allele Frequency") +
  theme_minimal()

# Arrange the two plots vertically
p1 + p2




library(dplyr)
library(ggplot2)
library(patchwork)

# Define the population size (e.g., births in 2023)
population_size <- 83702
# 
# # Calculate disease probability, expected cases, and probability of at least one case
# # for each variant, across all clinvar_clnsig categories.
# df_calc <- df %>%
#   # Convert allele frequency to numeric and remove missing values:
#   mutate(allele_freq = as.numeric(gnomAD_genomes_AF)) %>%
#   filter(!is.na(allele_freq)) %>%
#   # Compute disease probability based on Inheritance:
#   mutate(disease_prob = ifelse(Inheritance %in% c("AD", "X-linked"), 
#                                allele_freq,         # AD / X-linked: one copy is enough
#                                allele_freq^2),      # AR: need two copies
#          expected_cases = population_size * disease_prob,
#          prob_at_least_one = 1 - (1 - disease_prob)^population_size)


# Calculate disease probability, expected cases, and probability of at least one case
# for each variant, across all clinvar_clnsig categories.
df_calc <- df %>%
  mutate(allele_freq = as.numeric(gnomAD_genomes_AF)) %>%  # convert allele frequency to numeric
  filter(!is.na(allele_freq)) %>%
  mutate(disease_prob = ifelse(Inheritance %in% c("AD", "X-linked"), 
                               allele_freq,         # AD / X-linked: one copy is enough
                               allele_freq^2),      # AR: need two copies
         expected_cases = population_size * disease_prob,
         prob_at_least_one = 1 - (1 - disease_prob)^population_size)

# Define the full set of ClinVar categories (ensuring all labels are present)
clinvar_levels <- df$clinvar_clnsig |> unique()

# Tally the results per clinvar_clnsig category and complete missing ones with zeros
df_tally <- df_calc %>%
  group_by(clinvar_clnsig) %>%
  summarise(total_expected_cases = sum(expected_cases),
            overall_prob = 1 - prod(1 - disease_prob),
            .groups = "drop") %>%
  tidyr::complete(clinvar_clnsig = clinvar_levels,
           fill = list(total_expected_cases = 0, overall_prob = 0))

df_tally


# Plot 1: Expected cases vs allele frequency, with lines grouped by clinvar_clnsig
p1 <- ggplot(df_calc, aes(x = allele_freq, y = expected_cases, color = clinvar_clnsig)) +
  geom_point() +
  geom_line(aes(group = clinvar_clnsig)) +
  scale_x_log10() +
  labs(x = "Allele Frequency (log scale)",
       y = "Expected Cases",
       title = "Expected Cases by ClinVar Clinical Significance") +
  theme_minimal() +
  facet_grid(~clinvar_clnsig, scales = "free")
p1

# Plot 2: Probability of at least one case vs allele frequency, with lines grouped by clinvar_clnsig
p2 <- ggplot(df_calc, aes(x = allele_freq, y = prob_at_least_one, color = clinvar_clnsig)) +
  geom_point() +
  geom_line(aes(group = clinvar_clnsig)) +
  scale_x_log10() +
  labs(x = "Allele Frequency (log scale)",
       y = "Probability of ≥1 Case",
       title = "Probability of At Least One Case by ClinVar Clinical Significance") +
  theme_minimal() +
  facet_grid(~clinvar_clnsig, scales = "free")

p2

# Combine the plots using patchwork:
p1 / p2









# Total expected cases and overall probability across all variants:
total_expected_cases <- sum(df_calc$expected_cases)
overall_prob <- 1 - prod(1 - df_calc$disease_prob)

cat("Total expected cases (all variants):", total_expected_cases, "\n")
cat("Overall probability that a birth is affected by at least one variant:", overall_prob, "\n")
# 
# # Alternatively, get the tally per clinvar_clnsig category:
# df_tally <- df_calc %>%
#   group_by(clinvar_clnsig) %>%
#   summarise(total_expected_cases = sum(expected_cases),
#             overall_prob = 1 - prod(1 - disease_prob))
# 
# print(df_tally)

# Alternatively, get the tally per clinvar_clnsig category and ensure all levels are present:
clinvar_levels <- unique(df$clinvar_clnsig)

df_tally <- df_calc %>%
  group_by(clinvar_clnsig) %>%
  summarise(total_expected_cases = sum(expected_cases),
            overall_prob = 1 - prod(1 - disease_prob),
            .groups = "drop") %>%
  tidyr::complete(clinvar_clnsig = clinvar_levels,
                  fill = list(total_expected_cases = 0, overall_prob = 0))

print(df_tally)

# Plot the total expected cases by clinvar_clnsig as a bar chart, annotated with counts
p_bar <- ggplot(df_tally, aes(x = clinvar_clnsig, y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(total_expected_cases, 1)), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected Cases (per year)",
       title = "Total Expected Cases by\nClinVar Clinical Significance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Plot the overall probability (per category) as a bar chart, annotated with percentages
p_prob <- ggplot(df_tally, aes(x = clinvar_clnsig, y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.1)), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall Probability",
       title = "Overall Probability of an\nAffected Birth by ClinVar Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = scales::percent)

# Combine the two bar charts using patchwork:
p_bar + p_prob +
  plot_layout(guides = 'collect')
