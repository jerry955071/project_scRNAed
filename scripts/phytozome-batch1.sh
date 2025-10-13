#!/bin/bash

mkdir -p $1
cd $1
curl --cookie jgi_session=/api/sessions/a2dea57acdd5e6cfe3ad8263d7d6bd49 --output download.20250704.181929.zip -d "{\"ids\":{\"Phytozome-91\":{\"file_ids\":[\"52b9c8d1166e730e43a3503d\"],\"top_hit\":\"56901dbe0d878508e3d1fda2\"},\"Phytozome-297\":{\"file_ids\":[\"54ad8dde0d8785565d4707d3\"],\"top_hit\":\"54ad8dde0d8785565d4707d0\"},\"Phytozome-533\":{\"file_ids\":[\"5d94dc9ec0d65a87debccfcb\"],\"top_hit\":\"5d94dc9ec0d65a87debccfc8\"},\"Phytozome-870\":{\"file_ids\":[\"663e35f453447aa389b881a0\"],\"top_hit\":\"663e35f353447aa389b8818e\"}},\"api_version\":\"2\"}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
unzip *.zip