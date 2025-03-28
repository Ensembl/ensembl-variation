#!/bin/bash

####################
# Global Variables #
####################

# The Travis API endpoint. .com and .org are the commercial and free versions,
# respectively; enterprise users will have their own hostname.
endpoint=https://api.travis-ci.com

#############
# Functions #
#############

# Get this repo ID
repo_id () {
    curl -s -X GET -H "Authorization: token $AUTH_TOKEN" -H "Travis-API-Version: 3" https://api.travis-ci.com/repo/$1 | python3 -c "import sys, json; print(json.load(sys.stdin)['id'])"
}

# Make an API request using the auth token set above. First argument is the path
# of the API method, all later arguments are passed to curl directly.
travis_api () {
  curl -s $endpoint$1 \
       -H "Authorization: token $AUTH_TOKEN" \
       -H 'Content-Type: application/json' \
       -H 'Travis-API-Version: 3' \
       "${@:2}"
}

# Create a new environment variable for the repo and return its ID.
# First argument is the repo id, second is the environment variable
# name, and third is the value.
function env_var {
  travis_api /settings/env_vars?repository_id=$1 \
             -d "{\"env_var\":{\"name\":\"$2\",\"value\":\"$3\",\"public\":true}}" |
    sed 's/{"env_var":{"id":"\([^"]*\)",.*/\1/'
}

# print a spinner and terminate it
sp="/-\|"
sc=0
spin() {
   printf "\b${sp:sc++:1}"
   ((sc==${#sp})) && sc=0
}
endspin() {
   printf "\r%s\n" "$@"
}

# Only run for main builds. Pull request builds have the branch set to main,
# so ignore those too.
if [ "${TRAVIS_BRANCH}" != "main" ] || [ "${TRAVIS_PULL_REQUEST}" != "false" ]; then
  exit 0
fi

# The list of downstream dependent repos
dep_repos=("Ensembl%2Fensembl%2Dvep"
	   "Ensembl%2Fensembl%2Drest")

for dep_repo in "${dep_repos[@]}"; do
    # Get the ID of the dependent repo
    dep_repo_id=`repo_id $dep_repo`
    echo "Dependent repo: $dep_repo (ID: $dep_repo_id)"

    echo "Checking API triggered builds in the last hour"
    if travis_api /repo/$dep_repo/builds?build.event_type=api | python3 travisci/api_build_run_last_hour.py | grep -q "True"; then
	echo "Detected recent API-triggered build (run in the last hour) ... skip."
	continue
    fi
    
    echo "----------------------------------"
    echo "Triggering build on dependent repo"
    echo "----------------------------------"

    body="{
 \"request\": {
 \"message\": \"Build triggered by upstream $TRAVIS_REPO_SLUG repo (commit: $TRAVIS_COMMIT, branch: $TRAVIS_BRANCH).\",
 \"branch\": \"main\"
}}"

    # Make the request to trigger the build and get the ID of the request
    dep_repo_main_build_request_id=`travis_api /repo/$dep_repo/requests -H 'Accept: application/json' -X POST -d "$body" | python3 -c "import sys, json; print(json.load(sys.stdin)['request']['id'])"`
    echo "Build request ID: $dep_repo_main_build_request_id"

    # Wait until request is approved or max amount of time has passed
    i=0
    echo "Waiting for build request $dep_repo_main_build_request_id to be approved "
    build_request_approved=""
    until travis_api /repo/$dep_repo/request/$dep_repo_main_build_request_id | grep -q '"result": "approved"'; do
	spin
	sleep 5
	
	true $(( i++ ))
	if [ $i -eq 100 ]
	then
	    echo " reached max waiting time ... ABORT"
	    exit 1
	fi
    done
    endspin
    echo "Build request approved."

    # Get the ID of the main build.
    dep_repo_main_build_id=`travis_api /repo/$dep_repo/request/$dep_repo_main_build_request_id | python3 -c "import sys, json; print(json.load(sys.stdin)['builds'][0]['id'])"`
    echo "Build on $dep_repo main branch created (ID: $dep_repo_main_build_id)"

    # # Set the three environment variables needed, and capture their IDs so that they
    # # can be removed later.
    # env_var_ids=(`env_var $dep_repo_id DEPENDENT_BUILD true`
    #              `env_var $dep_repo_id TRIGGER_COMMIT $TRAVIS_COMMIT`
    #              `env_var $dep_repo_id TRIGGER_REPO $TRAVIS_REPO_SLUG`)
    
    # Wait for the build to start using the new environment variables.
    i=0
    printf "Waiting for build $dep_repo_main_build_id to start  "
    build_started=""
    until travis_api /build/$dep_repo_main_build_id | grep -q '"state": "started"'; do
	spin
	sleep 5
	
	true $(( i++ ))
	if [ $i -eq 100 ]
	then
	    echo " reached max waiting time ... stop waiting"
	    build_started="not yet"
	    break
	fi
    done
    endspin
    echo "Build $dep_repo_main_build_id $build_started started"

    # Remove all of the environment variables set above. This does mean that if this
    # script is terminated for whatever reason, these will need to be cleaned up
    # manually. We can do this either through the API, or by going to Settings ->
    # Environment Variables in the Travis web interface.
    #
    # for env_var_id in "${env_var_ids[@]}"; do
    #   travis_api /settings/env_vars/$env_var_id?repository_id=$dep_repo_id -X DELETE
    # done
done
