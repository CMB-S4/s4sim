#!/usr/bin/env python

# Based on
#   https://github.com/globus/automation-examples/blob/master/globus_folder_sync.py
# cloned on 01/05/2022

"""Sync a directory between two Globus endpoints.

Authorization only needs to happen once, afterwards tokens are saved to disk
(MUST BE STORED IN A SECURE LOCATION). Store data is already checked for
previous transfers, so if this script is run twice in quick succession,
the second run won't queue a duplicate transfer."""

import json
import sys
import os
import six

from globus_sdk import (NativeAppAuthClient, TransferClient,
                        RefreshTokenAuthorizer, TransferData,
                        GlobusAPIError, TransferAPIError)

from fair_research_login import NativeClient

# NERSC DTN endpoint
SOURCE_ENDPOINT = '9d6d994a-6d04-11e5-ba46-22000b92c6ec'

# NERSC DTN endpoint
#DESTINATION_ENDPOINT = '9d6d994a-6d04-11e5-ba46-22000b92c6ec'

# NERSC cmbs4 collaboration account endpoint
DESTINATION_ENDPOINT = 'fb077a06-6f06-11ec-b2c3-1b99bfd4976a'

# Copy data off of the endpoint share
#SOURCE_PATH = '/global/cscratch1/sd/keskital/s4sim/dc1/noise_sim/outputs/LAT0_CHLAT/f090/LAT0_CHLAT_split_schedule_1500'
SOURCE_PATH = '/global/cscratch1/sd/keskital/s4sim/dc1/cmb_sim/outputs'

# Destination Path -- The directory will be created if it doesn't exist
#DESTINATION_PATH = '/global/cfs/cdirs/cmbs4/dc/dc1/staging/noise_sim/outputs/LAT0_CHLAT/f090/LAT0_CHLAT_split_schedule_1500'
DESTINATION_PATH = '/global/cfs/cdirs/cmbs4/dc/dc1/staging/cmb_sim/outputs'

TRANSFER_LABEL = 'DC1 Sync Test'

# You will need to register a *Native App* at https://developers.globus.org/
# Your app should include the following:
#     - The scopes should match the SCOPES variable below
#     - Your app's clientid should match the CLIENT_ID var below
#     - "Native App" should be checked
# For more information:
# https://docs.globus.org/api/auth/developer-guide/#register-app
CLIENT_ID = '659aebe0-f505-46b5-8f13-348f451ef3d9'
DATA_FILE = 'transfer-data.json'
REDIRECT_URI = 'https://auth.globus.org/v2/web/auth-code'
SCOPES = ('openid email profile '
          'urn:globus:auth:scope:transfer.api.globus.org:all')

APP_NAME = 'globus_folder_sync.py'

# ONLY run new tasks if there was a previous task and it exited with one of the
# following statuses. This is ignored if there was no previous task.
# The previous task is queried from the DATA_FILE
PREVIOUS_TASK_RUN_CASES = ['SUCCEEDED', 'FAILED']

# Create the destination folder if it does not already exist
CREATE_DESTINATION_FOLDER = True


get_input = getattr(__builtins__, 'raw_input', input)


def load_data_from_file(filepath):
    """Load a set of saved tokens."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, 'r') as f:
        tokens = json.load(f)

    return tokens


def save_data_to_file(filepath, key, data):
    """Save data to a file"""
    try:
        store = load_data_from_file(filepath)
    except:
        store = {}
    if len(store) > 0:
        store[key] = data
    with open(filepath, 'w') as f:
        json.dump(store, f)


def setup_transfer_client(transfer_tokens):

    authorizer = RefreshTokenAuthorizer(
        transfer_tokens['refresh_token'],
        NativeAppAuthClient(client_id=CLIENT_ID),
        access_token=transfer_tokens['access_token'],
        expires_at=transfer_tokens['expires_at_seconds'])

    transfer_client = TransferClient(authorizer=authorizer)

    try:
        transfer_client.endpoint_autoactivate(SOURCE_ENDPOINT)
        transfer_client.endpoint_autoactivate(DESTINATION_ENDPOINT)
    except GlobusAPIError as ex:
        if ex.http_status == 401:
            sys.exit('Refresh token has expired. '
                     'Please delete the `tokens` object from '
                     '{} and try again.'.format(DATA_FILE))
        else:
            raise ex
    return transfer_client


def check_endpoint_path(transfer_client, endpoint, path):
    """Check the endpoint path exists"""
    try:
        transfer_client.operation_ls(endpoint, path=path)
    except TransferAPIError as tapie:
        print('Failed to query endpoint "{}": {}'.format(
            endpoint,
            tapie.message
        ))
        sys.exit(1)


def create_destination_directory(transfer_client, dest_ep, dest_path):
    """Create the destination path if it does not exist"""
    try:
        transfer_client.operation_ls(dest_ep, path=dest_path)
    except TransferAPIError:
        try:
            transfer_client.operation_mkdir(dest_ep, dest_path)
            print('Created directory: {}'.format(dest_path))
        except TransferAPIError as tapie:
            print('Failed to start transfer: {}'.format(tapie.message))
            sys.exit(1)


def main():
    tokens = None
    client = NativeClient(client_id=CLIENT_ID, app_name=APP_NAME)
    try:
        # if we already have tokens, load and use them
        tokens = client.load_tokens(requested_scopes=SCOPES)
    except:
        pass

    if not tokens:
        # if we need to get tokens, start the Native App authentication process
        # need to specify that we want refresh tokens
        tokens = client.login(requested_scopes=SCOPES,
                              refresh_tokens=True)
        try:
            client.save_tokens(tokens)
        except:
            pass

    transfer = setup_transfer_client(tokens['transfer.api.globus.org'])

    try:
        data = load_data_from_file(DATA_FILE)
        if len(data) > 0:
            task_data = data['task']
            task = transfer.get_task(task_data['task_id'])
            if task['status'] not in PREVIOUS_TASK_RUN_CASES:
                print('The last transfer status is {}, skipping run...'.format(
                    task['status']
                ))
                sys.exit(1)
    except KeyError:
        # Ignore if there is no previous task
        pass

    check_endpoint_path(transfer, SOURCE_ENDPOINT, SOURCE_PATH)
    if CREATE_DESTINATION_FOLDER:
        create_destination_directory(transfer, DESTINATION_ENDPOINT,
                                     DESTINATION_PATH)
    else:
        check_endpoint_path(transfer, DESTINATION_ENDPOINT, DESTINATION_PATH)

    tdata = TransferData(
        transfer,
        SOURCE_ENDPOINT,
        DESTINATION_ENDPOINT,
        label=TRANSFER_LABEL,
        sync_level="checksum"
    )
    tdata.add_item(SOURCE_PATH, DESTINATION_PATH, recursive=True)

    task = transfer.submit_transfer(tdata)
    save_data_to_file(DATA_FILE, 'task', task.data)
    print('Transfer has been started from\n  {}:{}\nto\n  {}:{}'.format(
        SOURCE_ENDPOINT,
        SOURCE_PATH,
        DESTINATION_ENDPOINT,
        DESTINATION_PATH
    ))
    url_string = 'https://app.globus.org/app/transfer?' + \
        six.moves.urllib.parse.urlencode({
            'origin_id': SOURCE_ENDPOINT,
            'origin_path': SOURCE_PATH,
            'destination_id': DESTINATION_ENDPOINT,
            'destination_path': DESTINATION_PATH
        })
    print('Visit the link below to see the changes:\n{}'.format(url_string))


if __name__ == '__main__':
    main()
