psql -U postgres -h localhost -c "drop schema zinc cascade;";
psql -U postgres -h localhost -c "delete from cgr_db_config where name='zinc'; create schema zinc;";
cgrdb create -c '{"host": "localhost", "password": "password", "user": "postgres"}' -n zinc -f config.json 

