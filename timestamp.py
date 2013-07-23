import datetime

def timestamp():
	return datetime.datetime.today().strftime("%a-%b-%d-%Y-%H:%M:%S")
