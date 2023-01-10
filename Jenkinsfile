#!/usr/bin/env groovy

//inspired by: https://www.christopherrung.com/2017/05/04/slack-build-notifications/

import groovy.json.JsonOutput

def slackNotificationChannel = "tp-commits"
def commit = ""
def author = ""
def modifiedFiles = ""

def notifySlack(text, channel, attachments) {

   def slackURL = sh(returnStdout: true, script: 'less /mn/stornext/d9/data/mikolajs/Jenkins_CI/webhook_bifrost.cf').trim()
   def payload = JsonOutput.toJson([text: text,
        channel: channel,
        username: "Jenkins",
        attachments: attachments
   ])
   sh "curl -X POST --data-urlencode \'payload=${payload}\' ${slackURL}"
}

def getGitAuthor = {
   commit = sh(returnStdout: true, script: 'git rev-parse --short HEAD')
   author = sh(returnStdout: true, script: "git --no-pager show -s --format='%an' ${commit}").trim()
   author = "`"+author+"`"
   commit = "`"+commit.trim()+"`"
}

def getModifiedFiles = {
   modifiedFiles = sh(returnStdout: true, script: 'git diff --name-only HEAD HEAD~1 | head -n 10').trim()
}

def uploadSummary(channel,file_name) {
   def botToken = sh(returnStdout: true, script: 'less /mn/stornext/d9/data/mikolajs/Jenkins_CI/app_bot_token.cf').trim()
   sh "curl -F file=@${file_name} -F \"initial_comment=Full report:\" -F channels=${channel} -H \"Authorization: Bearer ${botToken}\" https://slack.com/api/files.upload"
}

node {
   stage('Checkout'){
      checkout scm
   }

   try {

      stage("Post to Slack") {

         getGitAuthor()
         getModifiedFiles()

         // create output file with interesting information

         sh '''#!/bin/bash -el
         module load julia/1.8.5
	 export JULIA_LOAD_PATH=${JULIA_LOAD_PATH}:${PWD}/tests/modules
	 export JULIA_LOAD_PATH=${JULIA_LOAD_PATH}:${PWD}/src
	 julia tests/runUnitTests.jl | tee summary.out 
	 julia tests/experiments/testExBdrift.jl | tee -a summary.out
         '''

         def testSummary = sh(returnStdout: true, script:'cat summary.out').trim()

         testSummary = testSummary.replace("'","") //remove single quotes for now
         testSummary = "```" + testSummary + "```"

         notifySlack("", slackNotificationChannel, [
            [
               title: "${env.JOB_NAME}, build #${env.BUILD_NUMBER}",
               text: "Dear ${author}, congratulations your commit ${commit} passes my strict tests gracefully :+1:",
               color: "#00ff00",
               fields:[
                  [
                     title: "Modified Files",
                     value: "${modifiedFiles}",
                     short: true
                  ],
                  [
                     title: "Test Results",
                     value: "${testSummary}",
                     short: false
                  ]
               ]
            ]
         ])
      }

    } catch(e) {

      def testSummary = sh(returnStdout: true, script:'less summary.out').trim()

      testSummary = testSummary.replace("'","") //remove single quotes for now
      testSummary = "```" + testSummary + "```"

      notifySlack("", slackNotificationChannel, [
         [
            title: "${env.JOB_NAME}, build #${env.BUILD_NUMBER}",
            text: "Hello ${author} you might want to take a closer look to #${commit} commit.",
            color: "#ffd700",
            fields:[
               [
                  title: "Modified Files",
                  value: "${modifiedFiles}",
                  short: true
               ],
               [
                  title: "Test Results",
                  value: "${testSummary}",
                  short: false
               ]
            ]
         ]
      ])
      uploadSummary(slackNotificationChannel,'summary.std.txt')
   }
}


