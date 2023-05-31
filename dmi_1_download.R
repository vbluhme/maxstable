library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(geojsonsf)

dmi_stations <- GET(
  url = "https://dmigw.govcloud.dk/v2/climateData/collections/station/items?api-key=4848355b-8089-48a4-95b5-27d20a03688b"
) %>% 
  content("text") %>% 
  jsonlite::fromJSON(flatten = TRUE) %>% 
  .$features %>% 
  as_tibble() %>% 
  select(-type) %>% ## Otherwise naming duplicate with properties.type
  rename_with(~ str_remove(.x, "properties.")) %>% 
  mutate(
    x = map_dbl(geometry.coordinates, ~ .x[1]),
    y = map_dbl(geometry.coordinates, ~ .x[2])
  ) %>% 
  select(country, name, stationId, x, y) %>%
  group_by(stationId) %>% 
  filter(row_number() == n()) %>% 
  ungroup()

## Loads daily maximum temperatures for all weather stations for a given year.
load_dmi <- function(year) {
  url <- paste(
    "https://dmigw.govcloud.dk/v2/climateData/collections/stationValue/items?api-key=4848355b-8089-48a4-95b5-27d20a03688b",
    "parameterId=max_temp_w_date",
    "timeResolution=day",
    "limit=300000",
    paste0("datetime=", year, "-01-01T00:00:00Z/", year, "-12-31T00:00:00Z"),
    sep = "&"
  )
  
  GET(url = url) %>% 
    content("text") %>% 
    jsonlite::fromJSON(flatten = TRUE) %>% 
    .$features %>% 
    as_tibble() %>% 
    rename_with(~ str_remove(.x, "properties.")) %>% 
    rename(
      date = from,
      temp = value
    ) %>% 
    select(stationId, date, temp)
}

dmi <- tibble(year = 2011:2022) %>% 
  mutate(
    data = map(year, load_dmi)
  ) %>% 
  unnest(data) %>% 
  mutate(date = lubridate::as_date(date)) %>% 
  left_join(dmi_stations, by = "stationId") %>% 
  filter(
    country == "DNK",
    !(name %in% c("Nordby", "Kolding Lufthavn")), # Missing observations
    !(name %in% c("Bornholms Lufthavn", "Hammer Odde Fyr", "Nexø Vest")) # Bornholm
  ) %>% 
  select(stationId, year, date, temp)

dmi_stations <- dmi_stations %>% 
  filter(stationId %in% dmi$stationId) %>% 
  arrange(stationId)

save(dmi, dmi_stations, file = "dmi.RData")
