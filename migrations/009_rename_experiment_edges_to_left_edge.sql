-- Rename the experiment-level left-edge field to match its single value.
-- Stored values remain in raw-profile coordinates, without synthetic padding.
BEGIN IMMEDIATE;
ALTER TABLE experiments RENAME COLUMN EDGES TO LEFT_EDGE;
PRAGMA user_version=5;
COMMIT;
