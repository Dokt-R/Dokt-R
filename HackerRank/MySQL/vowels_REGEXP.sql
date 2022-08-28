-- Query the list of CITY names starting with vowels (i.e., a, e, i, o, or u) from STATION. Your result cannot contain duplicates.

-- SELECT DISTINCT city FROM station WHERE UPPER(LEFT(city,1)) IN ("A","E","I","O","U");
SELECT DISTINCT city FROM station WHERE city REGEXP '^[aeiou]'

-- Ending with vowels
SELECT DISTINCT city FROM station WHERE RIGHT(city,1) IN ("A","E","I","O","U")

-- Starting and ending
-- SELECT city FROM station WHERE UPPER(LEFT(city,1)) IN ("A","E","I","O","U") AND RIGHT(city,1) IN ("A","E","I","O","U")
SELECT DISTINCT city FROM station WHERE city REGEXP '^[aeiou].*[aeiou]$'

-- Not starting with vowels
SELECT DISTINCT city FROM station WHERE city NOT REGEXP '^[aeiou]'

-- Not ending
SELECT DISTINCT city FROM station WHERE city NOT REGEXP '[aeiou]$'

-- Not starting OR not ending
SELECT DISTINCT city FROM station WHERE city NOT REGEXP '^[aeiou].*[aeiou]$'

-- Not starting AND not ending
SELECT DISTINCT city FROM station WHERE city REGEXP '^[^aeiou].*.[^aeiou]$'