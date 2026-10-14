DROP TABLE IF EXISTS variants;
DROP TABLE IF EXISTS genes;

CREATE TABLE genes (
    gene TEXT PRIMARY KEY,
    role TEXT NOT NULL
);

CREATE TABLE variants (
    variation_id TEXT NOT NULL,
    assembly TEXT NOT NULL,
    gene TEXT NOT NULL,
    significance TEXT NOT NULL,
    review_stars TEXT NOT NULL,
    last_evaluated TEXT,
    start_pos INTEGER NOT NULL,
    stop_pos INTEGER NOT NULL,
    PRIMARY KEY (variation_id, assembly),
    CHECK (review_stars BETWEEN 0 AND 4),
    CHECK (stop_pos >= start_pos),
    FOREIGN KEY (gene) REFERENCES genes(gene)
);

CREATE INDEX idx_genes ON genes (gene);
CREATE INDEX idx_variants ON variants(significance);