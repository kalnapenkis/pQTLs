library(genpwr)
library(patchwork)
library(tidyverse)

pw_1 <- genpwr.calc(calc = "power",
                  model = "linear",
                  N = 500,
                  sd_y = 1,
                  k = NULL,
                  MAF = c(seq(0.001, 0.009, by = 0.001), seq(0.1, 0.5, by = 0.01)),
                  ES = 0.25,
                  Alpha = c(5e-8, 2.76e-10),
                  True.Model = c("Additive"),
                  Test.Model = c("Additive"))

pw_2 <- genpwr.calc(calc = "power",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    ES = 0.5,
                    Alpha = c(5e-8, 2.76e-10),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))


pw_3 <- genpwr.calc(calc = "power",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    ES = 1,
                    Alpha = c(5e-8, 2.76e-10),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))

pw_4 <- genpwr.calc(calc = "power",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    ES = 1.25,
                    Alpha = c(5e-8, 2.76e-10),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))


pw_5 <- genpwr.calc(calc = "power",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    ES = 1.5,
                    Alpha = c(5e-8, 2.76e-10),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))


pw <- rbind(pw_1, pw_2, pw_3, pw_4, pw_5)

colnames(pw)[5] <- "Effect size"

pw$`Effect size` <- as.character(pw$`Effect size`)

p1 <- ggplot(pw, aes(x = MAF, y = `Power_at_Alpha_5e-08`, colour = `Effect size`, group = `Effect size`)) +
  geom_line() + theme_classic() + scale_colour_brewer(palette = "Set2") +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "darkgreen") + ylab("Power") +
  ggtitle("P-value threshold 5e-8") + scale_x_continuous(breaks = c(seq(0, 0.09, by = 0.01), seq(0.1, 1, by = 0.1))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 <- ggplot(pw, aes(x = MAF, y = `Power_at_Alpha_2.76e-10`, colour = `Effect size`, group = `Effect size`)) +
  geom_line() + theme_classic() + scale_colour_brewer(palette = "Set2") +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "darkgreen") + ylab("Power") +
  ggtitle("P-value threshold 2.76e-10") + scale_x_continuous(breaks = c(seq(0, 0.09, by = 0.01), seq(0.1, 1, by = 0.1))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


p1 / p2

# calculate required effect size
sample_size_est_lenient <- genpwr.calc(calc = "es",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    Power = 0.8,
                    Alpha = c(5e-8),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))

sample_size_est_strict <- genpwr.calc(calc = "es",
                    model = "linear",
                    N = 500,
                    sd_y = 1,
                    k = NULL,
                    MAF = seq(0.001, 0.5, by = 0.01),
                    Power = 0.8,
                    Alpha = c(2.76e-10),
                    True.Model = c("Additive"),
                    Test.Model = c("Additive"))
