## Schema 28-29
## added new table, tagge_variation_feature
## that wiil contain variation_features that are
## tagged for a certain variation

CREATE TABLE tagged_variation_feature (

  variation_feature_id       INT not null,
  population_id              INT not null,
  
  PRIMARY KEY(variation_feature_id, population_id)
);
