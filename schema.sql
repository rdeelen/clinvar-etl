CREATE TABLE genes (
    gene TEXT PRIMARY KEY,
    role TEXT NOT NULL
);

CREATE TABLE variants (
    variation_id TEXT PRIMARY KEY,
    assembly TEXT PRIMARY KEY,
    significance TEXT NOT NULL,
    review_stars TEXT NOT NULL,
    last_evaluated TEXT NOT NULL,
    start_pos INTEGER NOT NULL,
    stop_pos INTEGER NOT NULL,
    CHECK (review_stars BETWEEN 0 AND 4),
    CHECK (stop_pos >= start_pos),
    FOREIGN KEY (gene)
        REFERENCES genes(gene)
);

CREATE INDEX index_genes
ON genes (gene)
ON variants(significance);