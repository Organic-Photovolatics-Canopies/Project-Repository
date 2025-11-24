`SHOW TABLES;`

```
+----------------------------+
| Tables_in_opv_data         |
+----------------------------+
| auth_group                 |
| auth_group_permissions     |
| auth_permission            |
| auth_user                  |
| auth_user_groups           |
| auth_user_user_permissions |
| data_calcqcset1            |
| data_calibqcset1           |
| data_geomscore             |
| data_geomstat              |
| data_graphscore            |
| data_graphstat             |
| data_molgeom               |
| data_molgraph              |
| data_scharber              |
| django_admin_log           |
| django_content_type        |
| django_session             |
| django_site                |
+----------------------------+
```

`DESCRIBE data_calcqcset1;`

```
+--------------------------+--------------+------+-----+---------+----------------+
| Field                    | Type         | Null | Key | Default | Extra          |
+--------------------------+--------------+------+-----+---------+----------------+
| id                       | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id             | int(11)      | NO   | MUL | NULL    |                |
| mol_geom_id              | int(11)      | NO   | MUL | NULL    |                |
| calc_id_str              | varchar(250) | NO   | UNI | NULL    |                |
| calc_tbz_str             | varchar(250) | NO   |     | NULL    |                |
| calc_archive_subdir_path | varchar(250) | NO   |     | NULL    |                |
| modelchem_str            | varchar(100) | NO   |     | NULL    |                |
| e_total                  | double       | YES  |     | NULL    |                |
| e_homo_alpha             | double       | YES  |     | NULL    |                |
| e_lumo_alpha             | double       | YES  |     | NULL    |                |
| e_gap_alpha              | double       | YES  |     | NULL    |                |
| e_homo_beta              | double       | YES  |     | NULL    |                |
| e_lumo_beta              | double       | YES  |     | NULL    |                |
| e_gap_beta               | double       | YES  |     | NULL    |                |
| e_gap_min                | double       | YES  |     | NULL    |                |
| dipmom_total             | double       | YES  |     | NULL    |                |
| s2_val                   | double       | YES  |     | NULL    |                |
+--------------------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_calibqcset1;`

```
+----------------+--------------+------+-----+---------+----------------+
| Field          | Type         | Null | Key | Default | Extra          |
+----------------+--------------+------+-----+---------+----------------+
| id             | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id   | int(11)      | YES  | MUL | NULL    |                |
| mol_geom_id    | int(11)      | YES  | MUL | NULL    |                |
| calc_qcset1_id | int(11)      | NO   | MUL | NULL    |                |
| calib_type     | varchar(100) | NO   |     | NULL    |                |
| e_homo_alpha   | double       | YES  |     | NULL    |                |
| e_lumo_alpha   | double       | YES  |     | NULL    |                |
| e_gap_alpha    | double       | YES  |     | NULL    |                |
| e_homo_beta    | double       | YES  |     | NULL    |                |
| e_lumo_beta    | double       | YES  |     | NULL    |                |
| e_gap_beta     | double       | YES  |     | NULL    |                |
| e_gap_min      | double       | YES  |     | NULL    |                |
+----------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_geomscore;`

```
+--------------+--------------+------+-----+---------+----------------+
| Field        | Type         | Null | Key | Default | Extra          |
+--------------+--------------+------+-----+---------+----------------+
| id           | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_geom_id  | int(11)      | NO   | MUL | NULL    |                |
| mol_graph_id | int(11)      | YES  | MUL | NULL    |                |
| score_type   | varchar(100) | NO   |     | NULL    |                |
| score        | double       | YES  |     | NULL    |                |
| score_n      | int(11)      | YES  |     | NULL    |                |
| score_min    | double       | YES  |     | NULL    |                |
| score_max    | double       | YES  |     | NULL    |                |
| score_mad    | double       | YES  |     | NULL    |                |
| score_rmsd   | double       | YES  |     | NULL    |                |
+--------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_geomstat;`

```
+---------------+--------------+------+-----+---------+----------------+
| Field         | Type         | Null | Key | Default | Extra          |
+---------------+--------------+------+-----+---------+----------------+
| id            | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id  | int(11)      | YES  | MUL | NULL    |                |
| mol_geom_id   | int(11)      | NO   | MUL | NULL    |                |
| calib_type    | varchar(100) | NO   |     | NULL    |                |
| property_type | varchar(100) | NO   |     | NULL    |                |
| average       | double       | YES  |     | NULL    |                |
| n             | int(11)      | YES  |     | NULL    |                |
| max           | double       | YES  |     | NULL    |                |
| min           | double       | YES  |     | NULL    |                |
| mad           | double       | YES  |     | NULL    |                |
| rmsd          | double       | YES  |     | NULL    |                |
+---------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_graphscore;`

```
+--------------+--------------+------+-----+---------+----------------+
| Field        | Type         | Null | Key | Default | Extra          |
+--------------+--------------+------+-----+---------+----------------+
| id           | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id | int(11)      | NO   | MUL | NULL    |                |
| score_type   | varchar(100) | NO   |     | NULL    |                |
| score        | double       | YES  |     | NULL    |                |
| score_n      | int(11)      | YES  |     | NULL    |                |
| score_min    | double       | YES  |     | NULL    |                |
| score_max    | double       | YES  |     | NULL    |                |
| score_mad    | double       | YES  |     | NULL    |                |
| score_rmsd   | double       | YES  |     | NULL    |                |
+--------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_graphstat;`

```
+---------------+--------------+------+-----+---------+----------------+
| Field         | Type         | Null | Key | Default | Extra          |
+---------------+--------------+------+-----+---------+----------------+
| id            | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id  | int(11)      | NO   | MUL | NULL    |                |
| calib_type    | varchar(100) | NO   |     | NULL    |                |
| property_type | varchar(100) | NO   |     | NULL    |                |
| average       | double       | YES  |     | NULL    |                |
| n             | int(11)      | YES  |     | NULL    |                |
| max           | double       | YES  |     | NULL    |                |
| min           | double       | YES  |     | NULL    |                |
| mad           | double       | YES  |     | NULL    |                |
| rmsd          | double       | YES  |     | NULL    |                |
+---------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_molgeom;`

```
+-------------------------+--------------+------+-----+---------+----------------+
| Field                   | Type         | Null | Key | Default | Extra          |
+-------------------------+--------------+------+-----+---------+----------------+
| id                      | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id            | int(11)      | NO   | MUL | NULL    |                |
| mol_id_str              | varchar(100) | NO   | UNI | NULL    |                |
| xyz_file_str            | varchar(100) | NO   |     | NULL    |                |
| xyz_archive_subdir_path | varchar(250) | NO   |     | NULL    |                |
| e_nucl                  | double       | YES  |     | NULL    |                |
| duplicate_geom          | tinyint(1)   | YES  |     | NULL    |                |
+-------------------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_molgraph;`

```
+------------------+--------------+------+-----+---------+----------------+
| Field            | Type         | Null | Key | Default | Extra          |
+------------------+--------------+------+-----+---------+----------------+
| id               | int(11)      | NO   | PRI | NULL    | auto_increment |
| SMILES_str       | varchar(250) | NO   | UNI | NULL    |                |
| iupac_str        | varchar(250) | NO   |     | NULL    |                |
| inchi_str        | varchar(250) | NO   |     | NULL    |                |
| cas_str          | varchar(250) | NO   |     | NULL    |                |
| trivial_namestr  | longtext     | NO   |     | NULL    |                |
| stoich_str       | varchar(30)  | NO   |     | NULL    |                |
| n_el             | int(11)      | YES  |     | NULL    |                |
| n_heavyatoms     | int(11)      | YES  |     | NULL    |                |
| n_bf_sz          | int(11)      | YES  |     | NULL    |                |
| n_bf_dzp         | int(11)      | YES  |     | NULL    |                |
| n_bf_tzp         | int(11)      | YES  |     | NULL    |                |
| mass             | double       | YES  |     | NULL    |                |
| permission_level | int(11)      | YES  |     | NULL    |                |
+------------------+--------------+------+-----+---------+----------------+
```

`DESCRIBE data_scharber;`

```
+---------------+--------------+------+-----+---------+----------------+
| Field         | Type         | Null | Key | Default | Extra          |
+---------------+--------------+------+-----+---------+----------------+
| id            | int(11)      | NO   | PRI | NULL    | auto_increment |
| mol_graph_id  | int(11)      | YES  | MUL | NULL    |                |
| mol_geom_id   | int(11)      | YES  | MUL | NULL    |                |
| calib_id      | int(11)      | YES  | MUL | NULL    |                |
| scharber_type | varchar(100) | NO   |     | NULL    |                |
| pce           | double       | YES  |     | NULL    |                |
| voc           | double       | YES  |     | NULL    |                |
| jsc           | double       | YES  |     | NULL    |                |
+---------------+--------------+------+-----+---------+----------------+
```
