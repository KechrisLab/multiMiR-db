CREATE USER 'multimir'@'localhost' IDENTIFIED BY 'multiMiR';
CREATE USER 'multimir'@'%' IDENTIFIED BY 'multiMiR';
GRANT ALL PRIVILEGES ON multimir.* TO 'multimir'@'localhost';
GRANT ALL PRIVILEGES ON multimir.* TO 'multimir'@'%';
