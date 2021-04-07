from . import util
from pathlib import Path


all_help = None # will be dict, once loaded

help_tags = ['table', 'figure']
help_directory = util.path_to_conga / 'help_messages'

def _load_help_messages():
    global all_help
    if all_help is None:
        all_help = {}
        for help_tag in help_tags:
            all_help[help_tag] = {}
            help_file = Path.joinpath(help_directory, help_tag+'s.txt')
            with open(help_file,'r') as data:
                current_tag = None
                for line in data:
                    if line.startswith('# ') and len(line.split()) == 2:
                        ## starting a new tag
                        current_tag = line.split()[1]
                        print('_load_help_messages: start new tag:',
                              current_tag)
                        all_help[help_tag][current_tag] = ''
                    elif current_tag is not None:
                        all_help[help_tag][current_tag] += line

# this is kind of silly
_load_help_messages()

def table_help(tag):
    if tag not in all_help['table']:
        return None
    return all_help['table'][tag]

def figure_help(tag):
    if tag not in all_help['figure']:
        return None
    return all_help['figure'][tag]

