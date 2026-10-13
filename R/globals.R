# Column symbols evaluated inside ggplot2/dplyr data masks, and the placeholder
# evaluated by scales::math_format(). Register only these intentional NSE names;
# ordinary functions are explicitly imported in NAMESPACE.
utils::globalVariables(c(
  ".x", "abundance", "classes", "conf.high", "conf.low", "effsize", "feature",
  "Freq", "group", "high", "importance", "low", "pd", "Prediction", "r", "rd",
  "Reference", "value", "variable", "x"
))
