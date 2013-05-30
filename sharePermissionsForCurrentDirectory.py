# sharePermissionsForCurrentDirectory.py
# Author: Brian Lin
# Description: Sets all files in the current directory to have read/write/execute permissions for user/group, and read permissions for world.

import os 

ls_output = os.listdir('.')
for entry in ls_output:
	if os.path.isdir(entry): continue
	os.chmod(entry, 00774)

