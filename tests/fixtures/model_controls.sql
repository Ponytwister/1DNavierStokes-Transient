CREATE TABLE model_controls (criterion TEXT PRIMARY KEY, value TEXT);
INSERT INTO model_controls VALUES
    ('width resolution (X)', '4'),
    ('length/time resolution (Z)', '2'),
    ('disable_reactions', 'true'),
    ('save_model_profiles', 'false'),
    ('max_iterations', '3'),
    ('convergence_epsx', '1e-9');
