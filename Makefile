all:
	pylint --disable=R0903,R0902,R1732 deapr.py
	pylint --disable=R0903,R0902,R1732,C0209 pathway.py
