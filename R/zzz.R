# TMB templates are compiled lazily in wnpmle_fit(), not during package load.
# This keeps library(wnpmle) lightweight and avoids runtime compilation in .onLoad().
