import datetime

def timestamp():
	return datetime.datetime.today().strftime("%a-%b-%d-%Y-%H-%M-%S")

def short():
	return datetime.datetime.today().strftime("%m.%d.%H%M")

def date():
	return datetime.datetime.today().strftime("%m.%d.%y")

def time():
	return datetime.datetime.today().strftime("%H-%M-%S")
